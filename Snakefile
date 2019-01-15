CONTRASTS = ['hd2o', 'd13acmw', 'd13d2o', 'd70acmw', 'd70d2o', 'd83acmw',
             'd83d2o']

FORCEFIELDS = ['martini', 'berger', 'slipids']
SURF_PRES = ['20', '30', '40', '50']

FIGURES_MAIN_TEXT = ['reports/figures/reflrefr.pdf',
                     'reports/figures/dspcdrywet.pdf',
                     'reports/figures/apm.pdf', 'reports/figures/trad_30.pdf',
                     'reports/figures/sim_slipids_30.pdf',
                     'reports/figures/sim_berger_30.pdf',
                     'reports/figures/sim_martini_30.pdf',
                     'reports/figures/number_density.pdf']

INPUTS_MAIN_TEXT = ['output/traditional/{}_30_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
INPUTS_MAIN_TEXT.append('output/simulation/{}_slipids_30_chisq.txt'.format(
    contrast) for contrast in CONTRASTS)
INPUTS_MAIN_TEXT.append('output/simulation/{}_berger_30_chisq.txt'.format(
    contrast) for contrast in CONTRASTS)
INPUTS_MAIN_TEXT.append('output/simulation/{}_martini_30_chisq.txt'.format(
    contrast) for contrast in CONTRASTS)
INPUTS_MAIN_TEXT.append('output/simulation/martini_30_tt.txt')
INPUTS_MAIN_TEXT.append('output/simulation/slipids_30_tt.txt')
INPUTS_MAIN_TEXT.append('output/simulation/berger_30_tt.txt')
INPUTS_MAIN_TEXT.append('output/simulation/slipids_30_wph.txt')
INPUTS_MAIN_TEXT.append('output/simulation/berger_30_wph.txt')
INPUTS_MAIN_TEXT.append('output/traditional/wph_30.txt')
INPUTS_MAIN_TEXT.append('output/traditional/ave_30_chisq.txt')

TRAD_REF_20 = ['output/traditional/{}_20_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
TRAD_SLD_20 = ['output/traditional/{}_20_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
TRAD_CHI_20 = ['output/traditional/{}_20_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
EXP_DATA_20 = ['data/experimental/surf_pres_20/{}20.dat'.format(
    contrast) for contrast in CONTRASTS]

TRAD_REF_30 = ['output/traditional/{}_30_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
TRAD_SLD_30 = ['output/traditional/{}_30_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
TRAD_CHI_30 = ['output/traditional/{}_30_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
EXP_DATA_30 = ['data/experimental/surf_pres_30/{}30.dat'.format(
    contrast) for contrast in CONTRASTS]

TRAD_REF_40 = ['output/traditional/{}_40_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
TRAD_SLD_40 = ['output/traditional/{}_40_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
TRAD_CHI_40 = ['output/traditional/{}_40_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
EXP_DATA_40 = ['data/experimental/surf_pres_40/{}40.dat'.format(
    contrast) for contrast in CONTRASTS]

TRAD_REF_50 = ['output/traditional/{}_50_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
TRAD_SLD_50 = ['output/traditional/{}_50_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
TRAD_CHI_50 = ['output/traditional/{}_50_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
EXP_DATA_50 = ['data/experimental/surf_pres_50/{}50.dat'.format(
    contrast) for contrast in CONTRASTS]

DENSITY_DATA = ['output/simulation/slipids_nb{}.txt'.format(
    index) for index in range(1, 11)]

SLI_DATA_20 = ['data/simulation/slipids/surf_pres_20/frame{}.pdb'.format(
    num) for num in range(1, 11)]
SLI_DATA_30 = ['data/simulation/slipids/surf_pres_30/frame{}.pdb'.format(
    num) for num in range(1, 11)]
SLI_DATA_40 = ['data/simulation/slipids/surf_pres_40/frame{}.pdb'.format(
    num) for num in range(1, 11)]
SLI_DATA_50 = ['data/simulation/slipids/surf_pres_50/frame{}.pdb'.format(
    num) for num in range(1, 11)]

BER_DATA_20 = ['data/simulation/berger/surf_pres_20/frame{}.pdb'.format(
    num) for num in range(1, 11)]
BER_DATA_30 = ['data/simulation/berger/surf_pres_30/frame{}.pdb'.format(
    num) for num in range(1, 11)]
BER_DATA_40 = ['data/simulation/berger/surf_pres_40/frame{}.pdb'.format(
    num) for num in range(1, 11)]
BER_DATA_50 = ['data/simulation/berger/surf_pres_50/frame{}.pdb'.format(
    num) for num in range(1, 11)]

MAR_DATA_20 = ['data/simulation/martini/surf_pres_20/frame{}.pdb'.format(
    num) for num in range(1, 11)]
MAR_DATA_30 = ['data/simulation/martini/surf_pres_30/frame{}.pdb'.format(
    num) for num in range(1, 11)]
MAR_DATA_40 = ['data/simulation/martini/surf_pres_40/frame{}.pdb'.format(
    num) for num in range(1, 11)]
MAR_DATA_50 = ['data/simulation/martini/surf_pres_50/frame{}.pdb'.format(
    num) for num in range(1, 11)]

SLI_REF_20 = ['output/simulation/{}_slipids_20_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_SLD_20 = ['output/simulation/{}_slipids_20_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_CHI_20 = ['output/simulation/{}_slipids_20_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_REF_30 = ['output/simulation/{}_slipids_30_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_SLD_30 = ['output/simulation/{}_slipids_30_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_CHI_30 = ['output/simulation/{}_slipids_30_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_REF_40 = ['output/simulation/{}_slipids_40_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_SLD_40 = ['output/simulation/{}_slipids_40_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_CHI_40 = ['output/simulation/{}_slipids_40_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_REF_50 = ['output/simulation/{}_slipids_50_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_SLD_50 = ['output/simulation/{}_slipids_50_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
SLI_CHI_50 = ['output/simulation/{}_slipids_50_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]

BER_REF_20 = ['output/simulation/{}_berger_20_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_SLD_20 = ['output/simulation/{}_berger_20_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_CHI_20 = ['output/simulation/{}_berger_20_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_REF_30 = ['output/simulation/{}_berger_30_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_SLD_30 = ['output/simulation/{}_berger_30_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_CHI_30 = ['output/simulation/{}_berger_30_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_REF_40 = ['output/simulation/{}_berger_40_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_SLD_40 = ['output/simulation/{}_berger_40_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_CHI_40 = ['output/simulation/{}_berger_40_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_REF_50 = ['output/simulation/{}_berger_50_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_SLD_50 = ['output/simulation/{}_berger_50_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
BER_CHI_50 = ['output/simulation/{}_berger_50_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]

MAR_REF_20 = ['output/simulation/{}_martini_20_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_SLD_20 = ['output/simulation/{}_martini_20_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_CHI_20 = ['output/simulation/{}_martini_20_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_REF_30 = ['output/simulation/{}_martini_30_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_SLD_30 = ['output/simulation/{}_martini_30_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_CHI_30 = ['output/simulation/{}_martini_30_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_REF_40 = ['output/simulation/{}_martini_40_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_SLD_40 = ['output/simulation/{}_martini_40_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_CHI_40 = ['output/simulation/{}_martini_40_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_REF_50 = ['output/simulation/{}_martini_50_ref.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_SLD_50 = ['output/simulation/{}_martini_50_sld.txt'.format(
    contrast) for contrast in CONTRASTS]
MAR_CHI_50 = ['output/simulation/{}_martini_50_chisq.txt'.format(
    contrast) for contrast in CONTRASTS]



rule all:
    input:
        'reports/preprint.pdf',
        'reports/si.pdf'

rule make_preprint:
    input:
        FIGURES_MAIN_TEXT,
        INPUTS_MAIN_TEXT,
        'reports/preprint.tex',
        'reports/paper.bib'
    output:
        'reports/preprint.pdf',
    shell:
        """
        cd reports && xelatex preprint.tex
        bibtex preprint.aux
        xelatex preprint.tex
        xelatex preprint.tex
        """

rule make_si:
    input:
        'reports/figures/trad_20.pdf',
        'reports/figures/trad_40.pdf',
        'reports/figures/trad_50.pdf',
        'reports/figures/sim_slipids_20.pdf',
        'reports/figures/sim_slipids_40.pdf',
        'reports/figures/sim_slipids_50.pdf',
        'output/simulation/ave_slipids_20_chisq.txt',
        'output/simulation/ave_slipids_40_chisq.txt',
        'output/simulation/ave_slipids_50_chisq.txt',
        'reports/figures/sim_berger_20.pdf',
        'reports/figures/sim_berger_40.pdf',
        'reports/figures/sim_berger_50.pdf',
        'output/simulation/ave_berger_20_chisq.txt',
        'output/simulation/ave_berger_40_chisq.txt',
        'output/simulation/ave_berger_50_chisq.txt',
        'output/simulation/slipids_30_tt.txt',
        'output/simulation/slipids_30_wph.txt',
        'output/simulation/martini_30_tt.txt',
        'output/simulation/martini_30_wph.txt',
        'reports/figures/sim_martini_20.pdf',
        'reports/figures/sim_martini_40.pdf',
        'reports/figures/sim_martini_50.pdf',
        'reports/figures/martiniorder.pdf',
        'reports/si.tex',
        'reports/paper.bib'
    output:
        'reports/si.pdf',
    shell:
        """
        cd reports && xelatex si.tex
        bibtex si.aux
        xelatex si.tex
        xelatex si.tex
        """

rule plot_apm:
    input:
        'bin/apm.py',
        'data/surf_iso.csv'
    output:
        'reports/figures/apm.pdf'
    shell:
        """
        cd bin && ipython apm.py
        cd ../
        """

rule trad_plot_20:
    input:
        TRAD_REF_20,
        TRAD_SLD_20,
        'notebooks/traditional/plot.py'
    output:
        'reports/figures/trad_20.pdf'
    run:
        sp = 20
        shell("cd notebooks/traditional && ipython plot.py {sp}")
        shell("cd ../")

rule trad_analysis_20:
    input:
        EXP_DATA_20,
        'notebooks/traditional/analysis.py',
        'models/mol_vol.py'
    output:
        TRAD_REF_20,
        TRAD_SLD_20,
        TRAD_CHI_20
    run:
        sp = 20
        shell("cd notebooks/traditional && ipython analysis.py {sp}")
        shell("cd ../")

rule trad_plot_30:
    input:
        TRAD_REF_30,
        TRAD_SLD_30,
        'notebooks/traditional/plot.py'
    output:
        'reports/figures/trad_30.pdf'
    run:
        sp = 30
        shell("cd notebooks/traditional && ipython plot.py {sp}")
        shell("cd ../")

rule trad_analysis_30:
    input:
        EXP_DATA_30,
        'notebooks/traditional/analysis.py',
        'models/mol_vol.py'
    output:
        TRAD_REF_30,
        TRAD_SLD_30,
        TRAD_CHI_30,
        'output/traditional/wph_30.txt'
    run:
        sp = 30
        shell("cd notebooks/traditional && ipython analysis.py {sp}")
        shell("cd ../")

rule trad_plot_40:
    input:
        TRAD_REF_40,
        TRAD_SLD_40,
        'notebooks/traditional/plot.py'
    output:
        'reports/figures/trad_40.pdf'
    run:
        sp = 40
        shell("cd notebooks/traditional && ipython plot.py {sp}")
        shell("cd ../")

rule trad_analysis_40:
    input:
        EXP_DATA_40,
        'notebooks/traditional/analysis.py',
        'models/mol_vol.py'
    output:
        TRAD_REF_40,
        TRAD_SLD_40,
        TRAD_CHI_40
    run:
        sp = 40
        shell("cd notebooks/traditional && ipython analysis.py {sp}")
        shell("cd ../")

rule trad_plot_50:
    input:
        TRAD_REF_50,
        TRAD_SLD_50,
        'notebooks/traditional/plot.py'
    output:
        'reports/figures/trad_50.pdf'
    run:
        sp = 50
        shell("cd notebooks/traditional && ipython plot.py {sp}")
        shell("cd ../")

rule trad_analysis_50:
    input:
        EXP_DATA_50,
        'notebooks/traditional/analysis.py',
        'models/mol_vol.py'
    output:
        TRAD_REF_50,
        TRAD_SLD_50,
        TRAD_CHI_50
    run:
        sp = 50
        shell("cd notebooks/traditional && ipython analysis.py {sp}")
        shell("cd ../")

rule trad_gen_analysis:
    input:
        'notebooks/traditional/analysis.ipynb'
    output:
        'notebooks/traditional/analysis.py'
    shell:
        """
        jupyter-nbconvert --to script {input}
        """

rule trad_gen_plot:
    input:
        'notebooks/traditional/plot.ipynb'
    output:
        'notebooks/traditional/plot.py'
    shell:
        """
        jupyter-nbconvert --to script {input}
        """

rule chisq_av:
    input:
        TRAD_CHI_20,
        TRAD_CHI_30,
        SLI_CHI_20,
        SLI_CHI_30,
        SLI_CHI_40,
        SLI_CHI_50,
        BER_CHI_20,
        BER_CHI_30,
        BER_CHI_40,
        BER_CHI_50,
        MAR_CHI_20,
        MAR_CHI_30,
        MAR_CHI_40,
        MAR_CHI_50,
        TRAD_CHI_40,
        TRAD_CHI_50,
        'bin/chisq_total.py'
    output:
        'output/traditional/ave_20_chisq.txt',
        'output/traditional/ave_30_chisq.txt',
        'output/traditional/ave_40_chisq.txt',
        'output/traditional/ave_50_chisq.txt',
        'output/simulation/ave_slipids_20_chisq.txt',
        'output/simulation/ave_slipids_30_chisq.txt',
        'output/simulation/ave_slipids_40_chisq.txt',
        'output/simulation/ave_slipids_50_chisq.txt',
        'output/simulation/ave_berger_20_chisq.txt',
        'output/simulation/ave_berger_30_chisq.txt',
        'output/simulation/ave_berger_40_chisq.txt',
        'output/simulation/ave_berger_50_chisq.txt'
    shell:
        """
        cd bin && ipython chisq_total.py
        cd ../
        """

rule nd_gen_plot:
    input:
        'notebooks/simulation/density_plot.ipynb'
    output:
        'notebooks/simulation/density_plot.py'
    shell:
        """
        jupyter-nbconvert --to script {input}
        """

rule nb_plot:
    input:
        DENSITY_DATA,
        'notebooks/simulation/density_plot.py'
    output:
        'reports/figures/number_density.pdf'
    run:
        shell("cd notebooks/simulation && ipython density_plot.py")
        shell("cd ../")

rule get_densities:
    input:
        SLI_DATA_30,
        'bin/get_density.py'
    output:
        DENSITY_DATA
    run:
        for index in range(1, 11):
            shell("cd bin && ipython get_density.py {index}")
            shell("cd ../")


rule sim_gen_analysis:
    input:
        'notebooks/simulation/analysis.ipynb'
    output:
        'notebooks/simulation/analysis.py'
    shell:
        """
        jupyter-nbconvert --to script {input}
        """

rule sim_gen_plot:
    input:
        'notebooks/simulation/plot.ipynb'
    output:
        'notebooks/simulation/plot.py'
    shell:
        """
        jupyter-nbconvert --to script {input}
        """

rule sli_plot_30:
    input:
        SLI_REF_30,
        SLI_SLD_30,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_slipids_30.pdf'
    run:
        ff = 'slipids'
        sp = '30'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule sli_analysis_30:
    input:
        EXP_DATA_30,
        SLI_DATA_30,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        SLI_REF_30,
        SLI_SLD_30,
        SLI_CHI_30
    run:
        ff = 'slipids'
        sp = '30'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")

rule sli_plot_20:
    input:
        SLI_REF_20,
        SLI_SLD_20,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_slipids_20.pdf'
    run:
        ff = 'slipids'
        sp = '20'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule sli_analysis_20:
    input:
        EXP_DATA_20,
        SLI_DATA_20,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        SLI_REF_20,
        SLI_SLD_20,
        SLI_CHI_20
    run:
        ff = 'slipids'
        sp = '20'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")

rule sli_plot_40:
    input:
        SLI_REF_40,
        SLI_SLD_40,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_slipids_40.pdf'
    run:
        ff = 'slipids'
        sp = '40'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule sli_analysis_40:
    input:
        EXP_DATA_40,
        SLI_DATA_40,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        SLI_REF_40,
        SLI_SLD_40,
        SLI_CHI_40
    run:
        ff = 'slipids'
        sp = '40'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")


rule sli_plot_50:
    input:
        SLI_REF_50,
        SLI_SLD_50,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_slipids_50.pdf'
    run:
        ff = 'slipids'
        sp = '50'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule sli_analysis_50:
    input:
        EXP_DATA_50,
        SLI_DATA_50,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        SLI_REF_50,
        SLI_SLD_50,
        SLI_CHI_50
    run:
        ff = 'slipids'
        sp = '50'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")

rule ber_plot_30:
    input:
        BER_REF_30,
        BER_SLD_30,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_berger_30.pdf'
    run:
        ff = 'berger'
        sp = '30'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule ber_analysis_30:
    input:
        EXP_DATA_30,
        BER_DATA_30,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        BER_REF_30,
        BER_SLD_30,
        BER_CHI_30
    run:
        ff = 'berger'
        sp = '30'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")

rule ber_plot_20:
    input:
        BER_REF_20,
        BER_SLD_20,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_berger_20.pdf'
    run:
        ff = 'berger'
        sp = '20'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule ber_analysis_20:
    input:
        EXP_DATA_20,
        BER_DATA_20,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        BER_REF_20,
        BER_SLD_20,
        BER_CHI_20
    run:
        ff = 'berger'
        sp = '20'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")

rule ber_plot_40:
    input:
        BER_REF_40,
        BER_SLD_40,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_berger_40.pdf'
    run:
        ff = 'berger'
        sp = '40'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule ber_analysis_40:
    input:
        EXP_DATA_40,
        BER_DATA_40,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        BER_REF_40,
        BER_SLD_40,
        BER_CHI_40
    run:
        ff = 'berger'
        sp = '40'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")


rule ber_plot_50:
    input:
        BER_REF_50,
        BER_SLD_50,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_berger_50.pdf'
    run:
        ff = 'berger'
        sp = '50'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule ber_analysis_50:
    input:
        EXP_DATA_50,
        BER_DATA_50,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        BER_REF_50,
        BER_SLD_50,
        BER_CHI_50
    run:
        ff = 'berger'
        sp = '50'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")

rule chain_tilt_gen:
    input:
        'notebooks/simulation/chain_tilt.ipynb'
    output:
        'notebooks/simulation/chain_tilt.py'
    shell:
        """
        jupyter-nbconvert --to script {input}
        """

rule martini_order_gen:
    input:
        'notebooks/simulation/martiniorder.ipynb'
    output:
        'notebooks/simulation/martiniorder.py'
    shell:
        """
        jupyter-nbconvert --to script {input}
        """

rule martini_order:
    input:
        'notebooks/simulation/martiniorder.py'
    output:
        'reports/figures/martiniorder.pdf'
    run:
        shell("cd notebooks/simulation && ipython martiniorder.py")


rule chain_tilt_berger:
    input:
        'notebooks/simulation/chain_tilt.py'
    output:
        'output/simulation/berger_30_tt.txt',
        'output/simulation/berger_30_wph.txt'
    run:
        ff = 'berger'
        for sp in SURF_PRES:
            sp = '_' + sp
            shell("cd notebooks/simulation && ipython chain_tilt.py {ff} {sp}")
            shell("cd ../")

rule chain_tilt_slipids:
    input:
        'notebooks/simulation/chain_tilt.py'
    output:
        'output/simulation/slipids_30_tt.txt',
        'output/simulation/slipids_30_wph.txt'
    run:
        ff = 'slipids'
        for sp in SURF_PRES:
            sp = '_' + sp
            shell("cd notebooks/simulation && ipython chain_tilt.py {ff} {sp}")
            shell("cd ../")

rule chain_tilt_martini:
    input:
        'notebooks/simulation/chain_tilt.py'
    output:
        'output/simulation/martini_30_tt.txt',
        'output/simulation/martini_30_wph.txt'
    run:
        ff = 'martini'
        for sp in SURF_PRES:
            sp = '_' + sp
            shell("cd notebooks/simulation && ipython chain_tilt.py {ff} {sp}")
            shell("cd ../")

rule mar_plot_30:
    input:
        MAR_REF_30,
        MAR_SLD_30,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_martini_30.pdf'
    run:
        ff = 'martini'
        sp = '30'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule mar_analysis_30:
    input:
        EXP_DATA_30,
        MAR_DATA_30,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        MAR_REF_30,
        MAR_SLD_30,
        MAR_CHI_30
    run:
        ff = 'martini'
        sp = '30'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")

rule mar_plot_20:
    input:
        MAR_REF_20,
        MAR_SLD_20,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_martini_20.pdf'
    run:
        ff = 'martini'
        sp = '20'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule mar_analysis_20:
    input:
        EXP_DATA_20,
        MAR_DATA_20,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        MAR_REF_20,
        MAR_SLD_20,
        MAR_CHI_20
    run:
        ff = 'martini'
        sp = '20'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")

rule mar_plot_40:
    input:
        MAR_REF_40,
        MAR_SLD_40,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_martini_40.pdf'
    run:
        ff = 'martini'
        sp = '40'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule mar_analysis_40:
    input:
        EXP_DATA_40,
        MAR_DATA_40,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        MAR_REF_40,
        MAR_SLD_40,
        MAR_CHI_40
    run:
        ff = 'martini'
        sp = '40'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")


rule mar_plot_50:
    input:
        MAR_REF_50,
        MAR_SLD_50,
        'notebooks/simulation/plot.py'
    output:
        'reports/figures/sim_martini_50.pdf'
    run:
        ff = 'martini'
        sp = '50'
        shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
        shell("cd ../")

rule mar_analysis_50:
    input:
        EXP_DATA_50,
        MAR_DATA_50,
        'notebooks/simulation/analysis.py',
        'bin/sim_lengths.py'
    output:
        MAR_REF_50,
        MAR_SLD_50,
        MAR_CHI_50
    run:
        ff = 'martini'
        sp = '50'
        for contrast in CONTRASTS:
            shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
            shell("cd ../")
