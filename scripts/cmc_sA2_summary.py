import pandas as pd
r = pd.read_csv('results_sim/cmc_sA2_sim.csv')
s = r.groupby('L').agg(reps=('rep', 'size'), n_sig_mean=('n_sig', 'mean'),
                       power=('power', 'mean'), fdr=('fdr', 'mean'),
                       fp_mean=('fp', 'mean')).round(3)
print('=== 按注入 LFC 汇总 ===')
print(s.to_string())
print('\n=== L=0 对照（应无显著）===')
print(r[r.L == 0][['rep', 'n_sig', 'fp']].to_string(index=False))
s.to_csv('results_sim/cmc_sA2_summary.csv')
