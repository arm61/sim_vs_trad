SURF_PRES = ['20', '30', '40', '50']
FORCEFIELDS = ['martini', 'berger', 'slipids']
CONTRASTS = ['d13acmw', 'd13d2o', 'hd2o', 'd70acmw', 'd70d2o', 'd83acmw', 'd83d2o']
METHODS = ['traditional', 'simulation', 'simulation', 'simulation']
FIGS_MAIN = ['reports/figures/reflrefr.pdf', 'reports/figures/dspcdrywet.png', 'reports/figures/apm.pdf', 'reports/figures/trad_30.pdf', 'reports/figures/sim_slipids_30.pdf', 'reports/figures/sim_berger_30.pdf', 'reports/figures/sim_martini_30.pdf', 'reports/figures/number_density.pdf']
ANAL_CHI = ['output/simulation/'+contrast+'_' + ff + '_'+sp+'_chisq.txt' for contrast in CONTRASTS for sp in SURF_PRES for ff in FORCEFIELDS]
FORCEFIELDS2 = ['', '_martini', '_berger', '_slipids']
TOTAL_CHI = ['output/'+m+'/ave'+FORCEFIELDS2[i]+'_'+sp+'_chisq.txt' for i, m in enumerate(METHODS) for sp in SURF_PRES]
OUTS_MAIN = ['output/simulation/martini_30_tt.txt', 'output/simulation/slipids_30_tt.txt', 'output/simulation/berger_30_tt.txt', 'output/simulation/slipids_30_wph.txt', 'output/simulation/berger_30_wph.txt', 'output/traditional/berger_30_wph.txt']

SURF_PRES2 = ['20', '40', '50']
FIGS_SI = ['reports/figures/martiniorder.pdf']
REFL_PLOTS_SI = ['reports/figures/sim_' + ff + '_' + sp + '.pdf' for sp in SURF_PRES2 for ff in FORCEFIELDS]
REFL_PLOTS_SI2 = ['reports/figures/trad_' + sp + '.pdf' for sp in SURF_PRES2]


#PARAMETERS = ['wph', 'dh', 'tt', 'angle']
#EXP_DATA = ['data/experimental/surf_pres_'+sp+'/'+contrast+sp+'.dat' for sp in SURF_PRES for #contrast in CONTRASTS]
#SIM_DATA = ['data/simulation/'+ff+'/surf_pres_'+sp+'/frame'+str(num)+'.pdb' for ff in FORCEFIELDS #for sp in SURF_PRES for num in range(1, 11)]
#CHAIN_TILT_OUT = ['output/simulation/'+ff+'_'+sp+'_'+para+'.txt' for ff in FORCEFIELDS for sp in #SURF_PRES for para in PARAMETERS]
#CHAIN_TILT_FIG = ['reports/figures/'+ff+'_'+sp+'_'+para+'.pdf' for ff in FORCEFIELDS for sp in #SURF_PRES for para in PARAMETERS]
#MAR_ANAL_REF = ['output/simulation/'+contrast+'_martini_'+sp+'_ref.txt' for contrast in CONTRASTS #for sp in SURF_PRES]
#MAR_ANAL_SLD = ['output/simulation/'+contrast+'_martini_'+sp+'_sld.txt' for contrast in CONTRASTS #for sp in SURF_PRES]
#BER_ANAL_REF = ['output/simulation/'+contrast+'_berger_'+sp+'_ref.txt' for contrast in CONTRASTS for #sp in SURF_PRES]
#BER_ANAL_SLD = ['output/simulation/'+contrast+'_berger_'+sp+'_sld.txt' for contrast in CONTRASTS for #sp in SURF_PRES]
#SLI_ANAL_REF = ['output/simulation/'+contrast+'_slipids_'+sp+'_ref.txt' for contrast in CONTRASTS #for sp in SURF_PRES]
#SLI_ANAL_SLD = ['output/simulation/'+contrast+'_slipids_'+sp+'_sld.txt' for contrast in CONTRASTS #for sp in SURF_PRES]
#SIM_FIGS = ['reports/figures/sim_'+ff+'_'+sp+'.pdf' for ff in FORCEFIELDS for sp in SURF_PRES]
#TRAD_ANAL_REF = ['output/traditional/'+contrast+'_'+sp+'_ref.txt' for contrast in CONTRASTS for sp #in SURF_PRES]
#TRAD_ANAL_SLD = ['output/traditional/'+contrast+'_'+sp+'_sld.txt' for contrast in CONTRASTS for sp #in SURF_PRES]
#TRAD_ANAL_CHI = ['output/traditional/'+contrast+'_'+sp+'_chisq.txt' for contrast in CONTRASTS for sp in SURF_PRES]#
#TRAD_FIGS = ['reports/figures/trad_'+sp+'.pdf' for sp in SURF_PRES]

#DENSITY_DATA = ['output/simulation/slipids_nb'+str(index)+'.txt' for index in range(1, 11)]

rule all:
    input:
        'reports/preprint.pdf',
        'reports/si.pdf'

rule make_preprint:
    input:
        'reports/figures/apm.pdf',
        #FIGS_MAIN,
        #ANAL_CHI,
        #TOTAL_CHI,
        #OUTS_MAIN,
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

rule makesi:
    input:
        #FIGS_SI,
        #REFL_PLOTS_SI,
        #REFL_PLOTS_SI2,
        'reports/si.tex',
        'reports/paper.bib'
    output:
        'reports/si.pdf'
    shell:
        """
        cd reports && xelatex si.tex
        bibtex si.aux
        xelatex si.tex
        xelatex si.tex
        cd ../
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

rule pdfclean:
    shell:
        "rm reports/paper.pdf"

rule clean:
    shell:
        "rm reports/paper.pdf"

"""rule sim_gen_analysis:
    input:
        'notebooks/simulation/analysis.ipynb'
    output:
        'notebooks/simulation/analysis.py'
    shell:

        jupyter-nbconvert --to script {input}


rule sim_gen_plot:
    input:
        'notebooks/simulation/plot.ipynb'
    output:
        'notebooks/simulation/plot.py'
    shell:

        jupyter-nbconvert --to script {input}


rule trad_gen_analysis:
    input:
        'notebooks/traditional/analysis.ipynb'
    output:
        'notebooks/traditional/analysis.py'
    shell:

        jupyter-nbconvert --to script {input}


rule trad_gen_analysis_mod:
    input:
        'notebooks/traditional/analysis_mod.ipynb'
    output:
        'notebooks/traditional/analysis_mod.py'
    shell:

        jupyter-nbconvert --to script {input}


rule trad_gen_plot:
    input:
        'notebooks/traditional/plot.ipynb'
    output:
        'notebooks/traditional/plot.py'
    shell:

        jupyter-nbconvert --to script {input}


rule nd_gen_plot:
    input:
        'notebooks/simulation/density_plot.ipynb'
    output:
        'notebooks/simulation/density_plot.py'
    shell:

        jupyter-nbconvert --to script {input}


rule chain_tilt_gen:
    input:
        'notebooks/simulation/chain_tilt.ipynb'
    output:
        'notebooks/simulation/chain_tilt.py'
    shell:

        jupyter-nbconvert --to script {input}


rule martini_analysis:
    input:
        EXP_DATA,
        SIM_DATA,
        'notebooks/simulation/analysis.py'
    output:
        MAR_ANAL_REF,
        MAR_ANAL_SLD,
        MAR_ANAL_CHI
    run:
        ff = 'martini'
        for sp in SURF_PRES:
            for contrast in CONTRASTS:
                shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
                shell("cd ../")

rule berger_analysis:
    input:
        EXP_DATA,
        SIM_DATA,
        'notebooks/simulation/analysis.py'
    output:
        BER_ANAL_REF,
        BER_ANAL_SLD,
        BER_ANAL_CHI
    run:
        ff = 'berger'
        for sp in SURF_PRES:
            for contrast in CONTRASTS:
                shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
                shell("cd ../")

rule slipids_analysis:
    input:
        EXP_DATA,
        SIM_DATA,
        'notebooks/simulation/analysis.py'
    output:
        SLI_ANAL_REF,
        SLI_ANAL_SLD,
        SLI_ANAL_CHI
    run:
        ff = 'slipids'
        for sp in SURF_PRES:
            for contrast in CONTRASTS:
                shell("cd notebooks/simulation && ipython analysis.py {ff} {sp} {contrast}")
                shell("cd ../")

rule sim_plot:
    input:
        MAR_ANAL_REF,
        MAR_ANAL_SLD,
        BER_ANAL_REF,
        BER_ANAL_SLD,
        SLI_ANAL_REF,
        SLI_ANAL_SLD,
        'notebooks/simulation/plot.py'
    output:
        SIM_FIGS
    run:
        for ff in FORCEFIELDS:
            for sp in SURF_PRES:
                shell("cd notebooks/simulation && ipython plot.py {ff} {sp}")
                shell("cd ../")

rule trad_analysis:
    input:
        EXP_DATA,
        'notebooks/traditional/analysis.py',
        'models/mol_vol.py'
    output:
        TRAD_ANAL_REF,
        TRAD_ANAL_SLD,
        TRAD_ANAL_CHI
    run:
        for sp in SURF_PRES:
            shell("cd notebooks/traditional && ipython analysis.py {sp}")
            shell("cd ../")

rule trad_plot:
    input:
        TRAD_ANAL_REF,
        TRAD_ANAL_SLD,
        'notebooks/traditional/plot.py'
    output:
        TRAD_FIGS
    run:
        for sp in SURF_PRES:
            shell("cd notebooks/traditional && ipython plot.py {sp}")
            shell("cd ../")

rule chisq_av:
    input:
        TRAD_ANAL_CHI,
        MAR_ANAL_CHI,
        BER_ANAL_CHI,
        SLI_ANAL_CHI,
        'bin/chisq_total.py'
    output:
        TOTAL_CHI
    shell:

        cd bin && ipython chisq_total.py
        cd ../


rule apm_plot:
    input:
        'bin/apm.py',
        'data/surf_iso.csv'
    output:
        'reports/figures/apm.pdf'
    shell:

        cd bin && ipython apm.py
        cd ../


rule get_densities:
    input:
        SIM_DATA,
        'bin/get_density.py'
    output:
        DENSITY_DATA
    run:
        for index in range(1, 11):
            shell("cd bin && ipython get_density.py {index}")
            shell("cd ../")

rule nb_plot:
    input:
        DENSITY_DATA,
        'notebooks/simulation/density_plot.py'
    output:
        'reports/figures/number_density.pdf'
    run:
        shell("cd notebooks/simulation && ipython density_plot.py")
        shell("cd ../")

rule chain_tilt:
    input:
        'notebooks/simulation/chain_tilt.py'
    output:
        CHAIN_TILT_FIG,
        CHAIN_TILT_OUT
    run:
        for ff in FORCEFIELDS:
            for sp in SURF_PRES:
                sp = '_' + sp
                shell("cd notebooks/simulation && ipython chain_tilt.py {ff} {sp}")
                shell("cd ../")"""
