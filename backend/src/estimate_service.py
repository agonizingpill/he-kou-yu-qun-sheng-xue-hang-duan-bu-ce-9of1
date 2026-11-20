import hashlib,json
from sqlalchemy import select
from .models import Transect,Echo,BioSample,LengthEstimate
from .common import finish,cruise_lock
from .school_service import latest
from . import length_model

def group_transects(transects):
    """航次航段按主航段整组归类：{根航段id:[航段,…]}。补航沿parent链归到主航段，已撤销航段不计。"""
    all_ts={t.id:t for t in transects}
    def root(t):
        while t.parent_id in all_ts:t=all_ts[t.parent_id]
        return t.id
    groups={}
    for t in transects:
        if t.status=='cancelled':continue
        groups.setdefault(root(t),[]).append(t)
    return groups

def basis_payload(analyses,samples):
    """采样依据的规范化表示：各航段最新已复核识别ID与对应生物样。补登或识别更新都会改变它。"""
    return {'analyses':{str(tid):a.id for tid,a in sorted(analyses.items())},'samples':[[s.id,s.analysis_id,s.cluster_label,s.count,s.mean_length_mm] for s in samples]}

def fingerprint(analyses,samples):
    return hashlib.sha256(json.dumps(basis_payload(analyses,samples),sort_keys=True,separators=(',',':')).encode('utf-8')).hexdigest()

def assemble(groups,analyses,samples,features):
    """纯数据组装估算结果。analyses:{航段id:最新已复核识别}；features:{识别id:{群落标签:{point_count,mean_strength_db,mean_depth_m}}}。
    返回 samples_used/group_errors/estimates/note，估算只覆盖未下网群落，不顶替实测。"""
    transects={t.id:t for members in groups.values() for t in members}
    group_of={t.id:root for root,members in groups.items() for t in members}
    by_analysis={a.id:(tid,a) for tid,a in analyses.items()}
    points,samples_used=[],[]
    for s in samples:
        tid,a=by_analysis[s.analysis_id]
        f=features.get(a.id,{}).get(s.cluster_label)
        if f is None:continue
        points.append({'key':s.id,'group':group_of[tid],'x':f['mean_strength_db'],'length_mm':s.mean_length_mm})
        samples_used.append({'sample_id':s.id,'analysis_id':a.id,'transect_id':tid,'transect':transects[tid].code,'group':transects[group_of[tid]].code,'cluster_label':s.cluster_label,'species':s.species,'count':s.count,'mean_length_mm':s.mean_length_mm,'point_count':f['point_count'],'mean_strength_db':round(f['mean_strength_db'],2),'mean_depth_m':round(f['mean_depth_m'],2)})
    targets,meta=[],[]
    for tid,a in sorted(analyses.items()):
        netted={s.cluster_label for s in samples if s.analysis_id==a.id}
        for label,f in sorted(features.get(a.id,{}).items()):
            if label in netted:continue
            targets.append({'key':len(targets),'group':group_of[tid],'x':f['mean_strength_db']})
            meta.append((tid,a.id,label,f))
    result=length_model.run(points,targets)
    group_errors=[{'group':transects[g['group']].code,'transects':[t.code for t in sorted(groups[g['group']],key=lambda t:t.id)],'sample_count':g['sample_count'],'mae_mm':round(g['mae_mm'],1),'clusters':[{'sample_id':c['key'],'actual_mm':round(c['actual_mm'],1),'predicted_mm':round(c['predicted_mm'],1),'error_mm':round(c['abs_error_mm'],1)} for c in g['clusters']]} for g in result['group_errors']]
    estimates=[]
    for e in result['estimates']:
        tid,aid,label,f=meta[e['key']]
        estimates.append({'transect_id':tid,'transect':transects[tid].code,'analysis_id':aid,'cluster_label':label,'point_count':f['point_count'],'mean_strength_db':round(f['mean_strength_db'],2),'mean_depth_m':round(f['mean_depth_m'],2),'mean_length_mm':round(e['value_mm'],1) if e['value_mm'] is not None else None,'blank_reason':e['blank_reason']})
    return samples_used,group_errors,estimates,result['note']

def _current_basis(db,cruise_id):
    """当前采样依据：整组分组、各航段最新已复核识别、这些识别对应的生物样。"""
    groups=group_transects(list(db.scalars(select(Transect).where(Transect.cruise_id==cruise_id))))
    analyses={}
    for members in groups.values():
        for t in members:
            a=latest(db,t.id,True)
            if a is not None:analyses[t.id]=a
    samples=list(db.scalars(select(BioSample).where(BioSample.analysis_id.in_([a.id for a in analyses.values()])).order_by(BioSample.id))) if analyses else []
    return groups,analyses,samples

def _features(db,analysis):
    """识别结果里各群落的回波点数、平均强度、平均水深。"""
    ids={}
    for p in analysis.assignments:
        if p['label']>=0:ids.setdefault(p['label'],[]).append(p['echo_id'])
    out={}
    for label,echo_ids in ids.items():
        rows=list(db.scalars(select(Echo).where(Echo.id.in_(echo_ids))))
        if rows:out[label]={'point_count':len(rows),'mean_strength_db':sum(e.strength_db for e in rows)/len(rows),'mean_depth_m':sum(e.depth_m for e in rows)/len(rows)}
    return out

def basis_sha256(db,cruise_id):
    _,analyses,samples=_current_basis(db,cruise_id)
    return fingerprint(analyses,samples)

def compute(db,cruise_id,user):
    c=cruise_lock(db,cruise_id)
    groups,analyses,samples=_current_basis(db,cruise_id)
    features={a.id:_features(db,a) for a in analyses.values()}
    samples_used,group_errors,estimates,note=assemble(groups,analyses,samples,features)
    row=LengthEstimate(cruise_id=c.id,basis_sha256=fingerprint(analyses,samples),note=note,samples_used=samples_used,group_errors=group_errors,estimates=estimates)
    db.add(row);db.flush();return finish(db,user,'估算未下网鱼群平均体长',row)
