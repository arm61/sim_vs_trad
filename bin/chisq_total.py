import numpy as np

def get_value(file, simtrad):
    f = open('../output/' + simtrad + '/' + file + '.txt', 'r')
    for line in f:
        k = line
    l = k.split('$')[1]
    return float(l)

simt = ['traditional', 'simulation', 'simulation', 'simulation']
ff = ['', '_martini', '_berger', '_slipids']
for j in range(0, len(simt)):
    sps = ['20', '30', '40', '50']
    for sp in sps:
        contrasts = ['d13acmw', 'd13d2o', 'd70acmw', 'd70d2o', 'd83acmw', 'd83d2o', 'hd2o']
        total = np.zeros(len(contrasts))
        for i in range(len(contrasts)):
            total[i] = get_value('{}{}_{}_chisq'.format(contrasts[i], ff[j], sp), simt[j])
        ave = total.mean()
        st = total.std()
        f_out = open('../output/{}/ave{}_{}_chisq.txt'.format(simt[j], ff[j], sp), 'w')
        f_out.write('${:.2f}\pm{:.2f}$'.format(ave, st))
