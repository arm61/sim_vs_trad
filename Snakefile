CONTRASTS = ['hd2o', 'd13acmw', 'd13d2o', 'd70acmw', 'd70d2o', 'd83acmw',
             'd83d2o']

FORCEFIELDS = ['martini', 'berger', 'slipids']
SURF_PRES = ['20', '30', '40', '50']

FIGURES_MAIN_TEXT = ['reports/figures/dspcdrywet.pdf',
                     'reports/figures/apm.pdf',
                     'reports/figures/dspc_30_ref_sld.pdf',
                     'reports/figures/dspc_slipids_30_ref_sld.pdf',
                     'reports/figures/dspc_berger_30_ref_sld.pdf',
                     'reports/figures/dspc_martini_30_ref_sld.pdf',
                     'reports/figures/dspc_slipids_30_ref_sld_short.pdf',
                     'reports/figures/water_30.pdf']

CHIS_MAIN_TEXT_TRAD = ['output/traditional/dspc_30_30_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_MAIN_TEXT_SLI = ['output/simulation/dspc_slipids_30_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_MAIN_TEXT_BER = ['output/simulation/dspc_berger_30_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_MAIN_TEXT_MAR = ['output/simulation/dspc_martini_30_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]

VARIABLES = ['-d_h', '-d_t', '_rough', '-phih', '-V_t', '-wph']

TRAD_MAIN_TEXT_20 = ['output/traditional/dspc_{}{}_{}.tex'.format('20', var, '20') for var in VARIABLES]
TRAD_MAIN_TEXT_30 = ['output/traditional/dspc_{}{}_{}.tex'.format('30', var, '30') for var in VARIABLES]
TRAD_MAIN_TEXT_40 = ['output/traditional/dspc_{}{}_{}.tex'.format('40', var, '40') for var in VARIABLES]
TRAD_MAIN_TEXT_50 = ['output/traditional/dspc_{}{}_{}.tex'.format('50', var, '50') for var in VARIABLES]

ATOMS = ['N', 'P', 'C2', 'C21', 'C31', 'C29', 'C39', '8C21', '8C31']
WHAT = ['mean', 'uq', 'position']

SPREAD_MAIN_TEXT = ['output/simulation/slipids_{}_{}_30.txt'.format(where, atom) for atom in ATOMS for where in WHAT]

INPUT_MAIN_TEXT = [CHIS_MAIN_TEXT_TRAD,
                   CHIS_MAIN_TEXT_SLI,
                   CHIS_MAIN_TEXT_BER,
                   CHIS_MAIN_TEXT_MAR,
                   TRAD_MAIN_TEXT_20,
                   TRAD_MAIN_TEXT_30,
                   TRAD_MAIN_TEXT_40,
                   TRAD_MAIN_TEXT_50,
                   'output/simulation/dspc_martini_30_dt.txt',
                   'output/simulation/dspc_slipids_30_dt.txt',
                   'output/simulation/dspc_berger_30_dt.txt',
                   'output/simulation/dspc_30_slipids_wph.txt',
                   'output/simulation/dspc_30_berger_wph.txt',
                   'output/traditional/dspc_30-wph_30.tex',
                   'output/simulation/slipids_position_C2_30.txt',
                   SPREAD_MAIN_TEXT]

FIGURES_SI_TEXT = ['reports/figures/martini_order.pdf',
                   'reports/figures/dspc_20_ref_sld.pdf',
                   'reports/figures/dspc_slipids_20_ref_sld.pdf',
                   'reports/figures/dspc_berger_20_ref_sld.pdf',
                   'reports/figures/dspc_martini_20_ref_sld.pdf',
                   'reports/figures/dspc_40_ref_sld.pdf',
                   'reports/figures/dspc_slipids_40_ref_sld.pdf',
                   'reports/figures/dspc_berger_40_ref_sld.pdf',
                   'reports/figures/dspc_martini_40_ref_sld.pdf',
                   'reports/figures/dspc_50_ref_sld.pdf',
                   'reports/figures/dspc_slipids_50_ref_sld.pdf',
                   'reports/figures/dspc_berger_50_ref_sld.pdf',
                   'reports/figures/dspc_martini_50_ref_sld.pdf',
                   'reports/figures/water_20.pdf',
                   'reports/figures/water_40.pdf',
                   'reports/figures/water_50.pdf']

SI_SURF_PRES = ['20', '30', '40', '50']

SPREAD_SI_TEXT_20 = ['output/simulation/slipids_{}_{}_{}.txt'.format(where, atom, '20') for atom in ATOMS for where in WHAT]
SPREAD_SI_TEXT_40 = ['output/simulation/slipids_{}_{}_{}.txt'.format(where, atom, '40') for atom in ATOMS for where in WHAT]
SPREAD_SI_TEXT_50 = ['output/simulation/slipids_{}_{}_{}.txt'.format(where, atom, '50') for atom in ATOMS for where in WHAT]

CHIS_SI_TEXT_TRAD_20 = ['output/traditional/dspc_20_20_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_SLI_20 = ['output/simulation/dspc_slipids_20_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_BER_20 = ['output/simulation/dspc_berger_20_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_MAR_20 = ['output/simulation/dspc_martini_20_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_TRAD_40 = ['output/traditional/dspc_40_40_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_SLI_40 = ['output/simulation/dspc_slipids_40_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_BER_40 = ['output/simulation/dspc_berger_40_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_MAR_40 = ['output/simulation/dspc_martini_40_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_TRAD_50 = ['output/traditional/dspc_50_50_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_SLI_50 = ['output/simulation/dspc_slipids_50_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_BER_50 = ['output/simulation/dspc_berger_50_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]
CHIS_SI_TEXT_MAR_50 = ['output/simulation/dspc_martini_50_{}_chi.txt'.format(
    contrast) for contrast in CONTRASTS]

INPUT_SI_TEXT = [SPREAD_SI_TEXT_20, SPREAD_SI_TEXT_40, SPREAD_SI_TEXT_50,
                 CHIS_SI_TEXT_TRAD_20, CHIS_SI_TEXT_SLI_20, CHIS_SI_TEXT_BER_20, CHIS_SI_TEXT_MAR_20,
                 CHIS_SI_TEXT_TRAD_40, CHIS_SI_TEXT_SLI_40, CHIS_SI_TEXT_BER_40, CHIS_SI_TEXT_MAR_40,
                 CHIS_SI_TEXT_TRAD_50, CHIS_SI_TEXT_SLI_50, CHIS_SI_TEXT_BER_50, CHIS_SI_TEXT_MAR_50]

EXP_DATA_20 = ['data/experimental/surf_pres_20/{}20.dat'.format(
    contrast) for contrast in CONTRASTS]
EXP_DATA_30 = ['data/experimental/surf_pres_30/{}30.dat'.format(
    contrast) for contrast in CONTRASTS]
EXP_DATA_40 = ['data/experimental/surf_pres_40/{}40.dat'.format(
    contrast) for contrast in CONTRASTS]
EXP_DATA_50 = ['data/experimental/surf_pres_50/{}50.dat'.format(
    contrast) for contrast in CONTRASTS]

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


rule all:
    input:
        'reports/preprint.pdf',
        'reports/si.pdf'

rule make_preprint:
    input:
        FIGURES_MAIN_TEXT,
        INPUT_MAIN_TEXT,
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
        FIGURES_SI_TEXT,
        INPUT_SI_TEXT,
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

rule sli_analysis_30:
    input:
        EXP_DATA_30,
        SLI_DATA_30,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_slipids_30_ref_sld.pdf',
        CHIS_MAIN_TEXT_SLI
    run:
        ff = 'slipids'
        sp = '30'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 1 0")
        shell("cd ../")

rule sli_wph_30:
    input:
        SLI_DATA_30,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_30_slipids_wph.txt'
    run:
        ff = 'slipids'
        sp = '30'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule sli_analysis_20:
    input:
        EXP_DATA_20,
        SLI_DATA_20,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_slipids_20_ref_sld.pdf',
        CHIS_SI_TEXT_SLI_20
    run:
        ff = 'slipids'
        sp = '20'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 1 0")
        shell("cd ../")

rule sli_wph_20:
    input:
        SLI_DATA_20,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_20_slipids_wph.txt'
    run:
        ff = 'slipids'
        sp = '20'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule sli_analysis_40:
    input:
        EXP_DATA_40,
        SLI_DATA_40,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_slipids_40_ref_sld.pdf',
        CHIS_SI_TEXT_SLI_40
    run:
        ff = 'slipids'
        sp = '40'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 1 0")
        shell("cd ../")

rule sli_wph_40:
    input:
        SLI_DATA_40,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_40_slipids_wph.txt'
    run:
        ff = 'slipids'
        sp = '40'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule sli_analysis_50:
    input:
        EXP_DATA_50,
        SLI_DATA_50,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_slipids_50_ref_sld.pdf',
        CHIS_SI_TEXT_SLI_50
    run:
        ff = 'slipids'
        sp = '50'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 1 0")
        shell("cd ../")

rule sli_wph_50:
    input:
        SLI_DATA_50,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_50_slipids_wph.txt'
    run:
        ff = 'slipids'
        sp = '50'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule ber_analysis_30:
    input:
        EXP_DATA_30,
        BER_DATA_30,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_berger_30_ref_sld.pdf',
        CHIS_MAIN_TEXT_BER
    run:
        ff = 'berger'
        sp = '30'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 1 0")
        shell("cd ../")

rule ber_wph_30:
    input:
        BER_DATA_30,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_30_berger_wph.txt'
    run:
        ff = 'berger'
        sp = '30'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule ber_analysis_20:
    input:
        EXP_DATA_20,
        BER_DATA_20,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_berger_20_ref_sld.pdf',
        CHIS_SI_TEXT_BER_20
    run:
        ff = 'berger'
        sp = '20'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 1 0")
        shell("cd ../")

rule ber_wph_20:
    input:
        BER_DATA_20,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_20_berger_wph.txt'
    run:
        ff = 'berger'
        sp = '20'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule ber_analysis_40:
    input:
        EXP_DATA_40,
        BER_DATA_40,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_berger_40_ref_sld.pdf',
        CHIS_SI_TEXT_BER_40
    run:
        ff = 'berger'
        sp = '40'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 1 0")
        shell("cd ../")

rule ber_wph_40:
    input:
        BER_DATA_40,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_40_berger_wph.txt'
    run:
        ff = 'berger'
        sp = '40'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule ber_analysis_50:
    input:
        EXP_DATA_50,
        BER_DATA_50,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_berger_50_ref_sld.pdf',
        CHIS_SI_TEXT_BER_50
    run:
        ff = 'berger'
        sp = '50'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 1 0")
        shell("cd ../")

rule ber_wph_50:
    input:
        BER_DATA_50,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_50_berger_wph.txt'
    run:
        ff = 'berger'
        sp = '50'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule mar_analysis_30:
    input:
        EXP_DATA_30,
        MAR_DATA_30,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_martini_30_ref_sld.pdf',
        CHIS_MAIN_TEXT_MAR
    run:
        ff = 'martini'
        sp = '30'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 4 0.4")
        shell("cd ../")

rule mar_wph_30:
    input:
        MAR_DATA_30,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_30_martini_wph.txt'
    run:
        ff = 'martini'
        sp = '30'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule mar_analysis_20:
    input:
        EXP_DATA_20,
        MAR_DATA_20,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_martini_20_ref_sld.pdf',
        CHIS_SI_TEXT_MAR_20
    run:
        ff = 'martini'
        sp = '20'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 4 0.4")
        shell("cd ../")

rule mar_wph_20:
    input:
        MAR_DATA_20,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_20_martini_wph.txt'
    run:
        ff = 'martini'
        sp = '20'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule mar_analysis_40:
    input:
        EXP_DATA_40,
        MAR_DATA_40,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_martini_40_ref_sld.pdf',
        CHIS_SI_TEXT_MAR_40
    run:
        ff = 'martini'
        sp = '40'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 4 0.4")
        shell("cd ../")

rule mar_wph_40:
    input:
        MAR_DATA_40,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_40_martini_wph.txt'
    run:
        ff = 'martini'
        sp = '40'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule mar_analysis_50:
    input:
        EXP_DATA_50,
        MAR_DATA_50,
        'scripts/simulation/md_simulation.py',
        'scripts/simulation/ref_help.py',
        'scripts/simulation/sim_lengths.py',
        'scripts/simulation/sim_analysis.py'
    output:
        'reports/figures/dspc_martini_50_ref_sld.pdf',
        CHIS_SI_TEXT_MAR_50
    run:
        ff = 'martini'
        sp = '50'
        shell("cd scripts/simulation && ipython sim_analysis.py {ff} {sp} 4 0.4")
        shell("cd ../")

rule mar_wph_50:
    input:
        MAR_DATA_50,
        'scripts/simulation/wph.py',
    output:
        'output/simulation/dspc_50_martini_wph.txt'
    run:
        ff = 'martini'
        sp = '50'
        shell("cd scripts/simulation && ipython wph.py {ff} {sp}")
        shell("cd ../")

rule trad_analysis_20:
    input:
        EXP_DATA_20,
        'scripts/traditional/chemically_consistent.py',
        'scripts/traditional/ref_help.py',
        'scripts/traditional/mol_vol.py'
    output:
        'output/traditional/dspc_20_chain.txt'
    run:
        shell("cd scripts/traditional && ipython chemically_consistent.py 20 47.9")
        shell("cd ../")

rule trad_analysis_30:
    input:
        EXP_DATA_30,
        'scripts/traditional/chemically_consistent.py',
        'scripts/traditional/ref_help.py',
        'scripts/traditional/mol_vol.py'
    output:
        'output/traditional/dspc_30_chain.txt'
    run:
        shell("cd scripts/traditional && ipython chemically_consistent.py 30 46.4")
        shell("cd ../")

rule trad_analysis_40:
    input:
        EXP_DATA_40,
        'scripts/traditional/chemically_consistent.py',
        'scripts/traditional/ref_help.py',
        'scripts/traditional/mol_vol.py'
    output:
        'output/traditional/dspc_40_chain.txt'
    run:
        shell("cd scripts/traditional && ipython chemically_consistent.py 40 45.0")
        shell("cd ../")

rule trad_analysis_50:
    input:
        EXP_DATA_50,
        'scripts/traditional/chemically_consistent.py',
        'scripts/traditional/ref_help.py',
        'scripts/traditional/mol_vol.py'
    output:
        'output/traditional/dspc_50_chain.txt'
    run:
        shell("cd scripts/traditional && ipython chemically_consistent.py 50 44.6")
        shell("cd ../")

rule trad_plot_20:
    input:
        'output/traditional/dspc_20_chain.txt',
        'scripts/traditional/cc_plot.py',
        'scripts/traditional/ref_help.py',
        'scripts/traditional/mol_vol.py'
    output:
        'reports/figures/dspc_20_pdf.pdf',
        'reports/figures/dspc_20_ref_sld.pdf',
        TRAD_MAIN_TEXT_20,
        CHIS_SI_TEXT_TRAD_20
    run:
        shell("cd scripts/traditional && ipython cc_plot.py 20 47.9")
        shell("cd ../")

rule trad_plot_30:
    input:
        'output/traditional/dspc_30_chain.txt',
        'scripts/traditional/cc_plot.py',
        'scripts/traditional/ref_help.py',
        'scripts/traditional/mol_vol.py'
    output:
        'reports/figures/dspc_30_pdf.pdf',
        'reports/figures/dspc_30_ref_sld.pdf',
        TRAD_MAIN_TEXT_30,
        CHIS_MAIN_TEXT_TRAD
    run:
        shell("cd scripts/traditional && ipython cc_plot.py 30 46.4")
        shell("cd ../")

rule trad_plot_40:
    input:
        'output/traditional/dspc_40_chain.txt',
        'scripts/traditional/cc_plot.py',
        'scripts/traditional/ref_help.py',
        'scripts/traditional/mol_vol.py'
    output:
        'reports/figures/dspc_40_pdf.pdf',
        'reports/figures/dspc_40_ref_sld.pdf',
        TRAD_MAIN_TEXT_40,
        CHIS_SI_TEXT_TRAD_40
    run:
        shell("cd scripts/traditional && ipython cc_plot.py 40 45.0")
        shell("cd ../")

rule trad_plot_50:
    input:
        'output/traditional/dspc_50_chain.txt',
        'scripts/traditional/cc_plot.py',
        'scripts/traditional/ref_help.py',
        'scripts/traditional/mol_vol.py'
    output:
        'reports/figures/dspc_50_pdf.pdf',
        'reports/figures/dspc_50_ref_sld.pdf',
        TRAD_MAIN_TEXT_50,
        CHIS_SI_TEXT_TRAD_50
    run:
        shell("cd scripts/traditional && ipython cc_plot.py 50 44.6")
        shell("cd ../")

rule waters_30:
    input:
        SLI_DATA_30,
        'scripts/simulation/waters.py'
    output:
        'output/simulation/waters_slipids_30.txt',
        'output/simulation/heads_slipids_30.txt',
        'output/simulation/tails_slipids_30.txt'
    run:
        shell("cd scripts/simulation && ipython waters.py slipids 30")
        shell("cd ../")

rule waters_plot_30:
    input:
        SLI_DATA_30,
        'output/simulation/waters_slipids_30.txt',
        'output/simulation/heads_slipids_30.txt',
        'output/simulation/tails_slipids_30.txt',
        'scripts/simulation/water_plot.py',
        'output/traditional/dspc_30-d_h_30.tex',
        'output/traditional/dspc_30-phih_30.tex'
    output:
        'reports/figures/water_30.pdf'
    run:
        shell("cd scripts/simulation && ipython water_plot.py 30")
        shell("cd ../")

rule waters_20:
    input:
        SLI_DATA_20,
        'scripts/simulation/waters.py'
    output:
        'output/simulation/waters_slipids_20.txt',
        'output/simulation/heads_slipids_20.txt',
        'output/simulation/tails_slipids_20.txt'
    run:
        shell("cd scripts/simulation && ipython waters.py slipids 20")
        shell("cd ../")

rule waters_plot_20:
    input:
        SLI_DATA_20,
        'output/simulation/waters_slipids_20.txt',
        'output/simulation/heads_slipids_20.txt',
        'output/simulation/tails_slipids_20.txt',
        'scripts/simulation/water_plot.py',
        'output/traditional/dspc_20-d_h_20.tex',
        'output/traditional/dspc_20-phih_20.tex'
    output:
        'reports/figures/water_20.pdf'
    run:
        shell("cd scripts/simulation && ipython water_plot.py 20")
        shell("cd ../")

rule waters_40:
    input:
        SLI_DATA_40,
        'scripts/simulation/waters.py'
    output:
        'output/simulation/waters_slipids_40.txt',
        'output/simulation/heads_slipids_40.txt',
        'output/simulation/tails_slipids_40.txt'
    run:
        shell("cd scripts/simulation && ipython waters.py slipids 40")
        shell("cd ../")

rule waters_plot_40:
    input:
        SLI_DATA_40,
        'output/simulation/waters_slipids_40.txt',
        'output/simulation/heads_slipids_40.txt',
        'output/simulation/tails_slipids_40.txt',
        'scripts/simulation/water_plot.py',
        'output/traditional/dspc_40-d_h_40.tex',
        'output/traditional/dspc_40-phih_40.tex'
    output:
        'reports/figures/water_40.pdf'
    run:
        shell("cd scripts/simulation && ipython water_plot.py 40")
        shell("cd ../")

rule waters_50:
    input:
        SLI_DATA_50,
        'scripts/simulation/waters.py'
    output:
        'output/simulation/waters_slipids_50.txt',
        'output/simulation/heads_slipids_50.txt',
        'output/simulation/tails_slipids_50.txt'
    run:
        shell("cd scripts/simulation && ipython waters.py slipids 50")
        shell("cd ../")

rule waters_plot_50:
    input:
        SLI_DATA_50,
        'output/simulation/waters_slipids_50.txt',
        'output/simulation/heads_slipids_50.txt',
        'output/simulation/tails_slipids_50.txt',
        'scripts/simulation/water_plot.py',
        'output/traditional/dspc_50-d_h_50.tex',
        'output/traditional/dspc_50-phih_50.tex'
    output:
        'reports/figures/water_50.pdf'
    run:
        shell("cd scripts/simulation && ipython water_plot.py 50")
        shell("cd ../")

rule spread_30:
    input:
        SLI_DATA_30,
        'scripts/simulation/roughness.py'
    output:
        SPREAD_MAIN_TEXT
    run:
        shell("cd scripts/simulation && ipython roughness.py 30")
        shell("cd ../../")


rule spread_20:
    input:
        SLI_DATA_20,
        'scripts/simulation/roughness.py'
    output:
        SPREAD_SI_TEXT_20
    run:
        shell("cd scripts/simulation && ipython roughness.py 20")
        shell("cd ../../")


rule spread_40:
    input:
        SLI_DATA_40,
        'scripts/simulation/roughness.py'
    output:
        SPREAD_SI_TEXT_40
    run:
        shell("cd scripts/simulation && ipython roughness.py 40")
        shell("cd ../../")


rule spread_50:
    input:
        SLI_DATA_50,
        'scripts/simulation/roughness.py'
    output:
        SPREAD_SI_TEXT_50
    run:
        shell("cd scripts/simulation && ipython roughness.py 50")
        shell("cd ../../")

rule martiniorder:
    input:
        MAR_DATA_30,
        'scripts/simulation/martini_order.py'
    output:
        'reports/figures/martini_order.pdf'
    run:
        shell("cd scripts/simulation && ipython martini_order.py")

rule tail_length_sli:
    input:
        SLI_DATA_30,
        'scripts/simulation/tail_length.py'
    output:
        'output/simulation/dspc_slipids_30_dt.txt'
    run:
        shell("cd scripts/simulation && ipython tail_length.py slipids 30")

rule tail_length_ber:
    input:
        BER_DATA_30,
        'scripts/simulation/tail_length.py'
    output:
        'output/simulation/dspc_berger_30_dt.txt'
    run:
        shell("cd scripts/simulation && ipython tail_length.py berger 30")

rule tail_length_mar:
    input:
        MAR_DATA_30,
        'scripts/simulation/tail_length.py'
    output:
        'output/simulation/dspc_martini_30_dt.txt'
    run:
        shell("cd scripts/simulation && ipython tail_length.py martini 30")

rule short_sli:
    input:
        SLI_DATA_30,
        'scripts/simulation/short_sim_analysis.py'
    output:
        'reports/figures/dspc_slipids_30_ref_sld_short.pdf'
    run:
        ff = 'slipids'
        sp = '30'
        shell("cd scripts/simulation && ipython short_sim_analysis.py {ff} {sp} 1 0")
        shell("cd ../")

rule clean:
    shell:
        """
        rm reports/si.pdf
        rm reports/preprint.pdf
        rm reports/figures/dspc_*.pdf
        rm reports/figures/water_*.pdf
        rm reports/figures/martini_order.pdf
        rm output/simulation/*
        rm output/traditional/*
        """
