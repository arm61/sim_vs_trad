SURF_PRES = ['20', '30', '40', '50']
CONTRASTS = ['d13acmw', 'd13d2o', 'hd2o', 'd70acmw', 'd70d2o', 'd83acmw', 'd83d2o']
EXP_DATA = ['data/experimental/surf_pres_'+sp+'/'+contrast+sp+'.dat' for sp in SURF_PRES for contrast in CONTRASTS]
FORCEFIELDS = ['martini', 'berger', 'slipids']
SIM_DATA = ['data/simulation/'+ff+'/surf_pres_'+sp+'/frame'+str(num)+'.pdb' for ff in FORCEFIELDS for sp in SURF_PRES for num in range(1, 11)]
MAR_ANAL_REF = ['output/simulation/'+contrast+'_martini_'+sp+'_ref.txt' for contrast in CONTRASTS for sp in SURF_PRES]
MAR_ANAL_SLD = ['output/simulation/'+contrast+'_martini_'+sp+'_sld.txt' for contrast in CONTRASTS for sp in SURF_PRES]
MAR_ANAL_CHI = ['output/simulation/'+contrast+'_martini_'+sp+'_chisq.txt' for contrast in CONTRASTS for sp in SURF_PRES]
BER_ANAL_REF = ['output/simulation/'+contrast+'_berger_'+sp+'_ref.txt' for contrast in CONTRASTS for sp in SURF_PRES]
BER_ANAL_SLD = ['output/simulation/'+contrast+'_berger_'+sp+'_sld.txt' for contrast in CONTRASTS for sp in SURF_PRES]
BER_ANAL_CHI = ['output/simulation/'+contrast+'_berger_'+sp+'_chisq.txt' for contrast in CONTRASTS for sp in SURF_PRES]
SLI_ANAL_REF = ['output/simulation/'+contrast+'_slipids_'+sp+'_ref.txt' for contrast in CONTRASTS for sp in SURF_PRES]
SLI_ANAL_SLD = ['output/simulation/'+contrast+'_slipids_'+sp+'_sld.txt' for contrast in CONTRASTS for sp in SURF_PRES]
SLI_ANAL_CHI = ['output/simulation/'+contrast+'_slipids_'+sp+'_chisq.txt' for contrast in CONTRASTS for sp in SURF_PRES]
SIM_FIGS = ['reports/figures/sim_'+ff+'_'+sp+'.pdf' for ff in FORCEFIELDS for sp in SURF_PRES]
TRAD_ANAL_REF = ['output/traditional/'+contrast+'_'+sp+'_ref.txt' for contrast in CONTRASTS for sp in SURF_PRES]
#MOD_ANAL_REF = ['output/traditional/'+contrast+'_'+sp+'_mod_ref.txt' for contrast in CONTRASTS for sp in SURF_PRES]
TRAD_ANAL_SLD = ['output/traditional/'+contrast+'_'+sp+'_sld.txt' for contrast in CONTRASTS for sp in SURF_PRES]
#MOD_ANAL_SLD = ['output/traditional/'+contrast+'_'+sp+'_mod_sld.txt' for contrast in CONTRASTS for sp in SURF_PRES]
TRAD_ANAL_CHI = ['output/traditional/'+contrast+'_'+sp+'_chisq.txt' for contrast in CONTRASTS for sp in SURF_PRES]
#MOD_ANAL_CHI = ['output/traditional/'+contrast+'_'+sp+'_mod_chisq.txt' for contrast in CONTRASTS for sp in SURF_PRES]
TRAD_FIGS = ['reports/figures/trad_'+sp+'.pdf' for sp in SURF_PRES]
#MOD_FIGS = ['reports/figures/trad_'+sp+'_mod.pdf' for sp in SURF_PRES]
METHODS = ['traditional', 'simulation', 'simulation', 'simulation']
FORCEFIELDS2 = ['', '_martini', '_berger', '_slipids']
TOTAL_CHI = ['output/'+m+'/ave'+FORCEFIELDS2[i]+'_'+sp+'_chisq.txt' for i, m in enumerate(METHODS) for sp in SURF_PRES]
DENSITY_DATA = ['output/simulation/slipids_nb'+str(index)+'.txt' for index in range(1, 11)]

rule all:
    input:
        'reports/paper.pdf'

rule make_pdf:
    input:
        SIM_FIGS,
        TRAD_FIGS,
#        MOD_FIGS,
        TOTAL_CHI,
        'reports/figures/apm.pdf',
        'reports/figures/number_density.pdf',
        'reports/paper.tex',
        'reports/paper.bib'
    output:
        'reports/paper.pdf'
    shell:
        """
        cd reports && xelatex paper.tex
        bibtex paper.aux
        xelatex paper.tex
        xelatex paper.tex
        cd ../
        """

rule pdfclean:
    shell:
        "rm reports/paper.pdf"

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

rule trad_gen_analysis:
    input:
        'notebooks/traditional/analysis.ipynb'
    output:
        'notebooks/traditional/analysis.py'
    shell:
        """
        jupyter-nbconvert --to script {input}
        """

rule trad_gen_analysis_mod:
    input:
        'notebooks/traditional/analysis_mod.ipynb'
    output:
        'notebooks/traditional/analysis_mod.py'
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

rule nd_gen_plot:
    input:
        'notebooks/simulation/density_plot.ipynb'
    output:
        'notebooks/simulation/density_plot.py'
    shell:
        """
        jupyter-nbconvert --to script {input}
        """

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

#rule trad_analysis_mod:
#    input:
#        EXP_DATA,
#        'notebooks/traditional/analysis_mod.py',
#        'models/mol_vol_mod.py'
#    output:
#        MOD_ANAL_REF,
#        MOD_ANAL_SLD,
#        MOD_ANAL_CHI
#    run:
#        for sp in SURF_PRES:
#            shell("cd notebooks/traditional && ipython analysis_mod.py {sp}")
#            shell("cd ../")
#
#rule trad_plot_mod:
#    input:
#        MOD_ANAL_REF,
#        MOD_ANAL_SLD,
#        'notebooks/traditional/plot.py'
#    output:
#        MOD_FIGS
#    run:
#        sp = ['20_mod', '30_mod', '40_mod', '50_mod']
#        for i in sp:
#            shell("cd notebooks/traditional && ipython plot.py {i}")
#            shell("cd ../")

rule chisq_av:
    input:
        TRAD_ANAL_CHI,
        MAR_ANAL_CHI,
        BER_ANAL_CHI,
        SLI_ANAL_CHI,
#        MOD_ANAL_CHI,
        'bin/chisq_total.py'
    output:
        TOTAL_CHI
    shell:
        """
        cd bin && ipython chisq_total.py
        cd ../
        """

rule apm_plot:
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
