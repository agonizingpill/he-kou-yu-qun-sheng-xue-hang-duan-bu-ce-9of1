"""未下网鱼群平均体长估算的纯计算核心。

声学依据：鱼体目标强度（dB）与体长对数近似线性，故用下网样本群落的
回波平均强度对 log10(平均体长mm) 做单变量最小二乘拟合，再反推未下网群落。
留验以"主航段及其补航"整组为单位留出，逐组验证并记录毫米误差。
只用标准库，便于脱离数据库独立复核。
"""
import math


def fit(xs, ys):
    """最小二乘 y=a+b*x；自变量无变异（或单点）时退化为常数模型。"""
    n = len(xs)
    if n == 0:
        raise ValueError('拟合至少需要一个样本点')
    if len(set(xs)) < 2:
        return (sum(ys) / n, 0.0)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    return (my - slope * mx, slope)


def _predict(model, x):
    a, b = model
    return a + b * x


def run(samples, targets, min_groups=2):
    """samples: [{'key','group','x','length_mm'}] 下网样本点，x 为群落回波平均强度 dB。
    targets: [{'key','group','x'}] 未下网群落，estimates 按相同顺序返回。
    返回 {'note','group_errors','estimates'}；留验误差与估算值均为毫米。
    可整组留验的组不足 min_groups 时全部留空；目标强度超出样本范围时该目标留空。"""
    points = [dict(p, y=math.log10(p['length_mm'])) for p in samples]
    groups = sorted({p['group'] for p in points})
    result = {'note': '', 'group_errors': [], 'estimates': []}
    if len(groups) < min_groups:
        result['note'] = '可整组留验的航段组不足两组，体长估算留空'
        result['estimates'] = [{'key': t['key'], 'value_mm': None, 'blank_reason': result['note']} for t in targets]
        return result
    for group in groups:
        train = [p for p in points if p['group'] != group]
        model = fit([p['x'] for p in train], [p['y'] for p in train])
        clusters = []
        for p in [p for p in points if p['group'] == group]:
            predicted = 10 ** _predict(model, p['x'])
            clusters.append({'key': p['key'], 'actual_mm': p['length_mm'], 'predicted_mm': predicted, 'abs_error_mm': abs(predicted - p['length_mm'])})
        result['group_errors'].append({'group': group, 'sample_count': len(clusters), 'mae_mm': sum(c['abs_error_mm'] for c in clusters) / len(clusters), 'clusters': clusters})
    model = fit([p['x'] for p in points], [p['y'] for p in points])
    lo = min(p['x'] for p in points)
    hi = max(p['x'] for p in points)
    for t in targets:
        if t['x'] < lo or t['x'] > hi:
            result['estimates'].append({'key': t['key'], 'value_mm': None, 'blank_reason': '目标回波强度超出样本范围，体长估算留空'})
        else:
            result['estimates'].append({'key': t['key'], 'value_mm': 10 ** _predict(model, t['x']), 'blank_reason': ''})
    return result
