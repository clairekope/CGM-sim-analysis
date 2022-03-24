# Scripts are run in paper figure order

# dependencies: ../sample_data/fid     DD0000,DD0020,DD0040,DD0060,DD0080
#               ../sample_data/cflow   DD0060
#               ../sample_data/tctff5  DD0060
#               ../sample_data/tctff20 DD0060
#               ../sample_data/linrot  DD0060
#               ../sample_data/norot   DD0060
#               ../extracted_data/{sim}_masses_over_time.txt
#               ../extracted_data/{sim}_sfh.txt
#               ../extracted_data/{sim}_tcool_mass_dist_CGM.txt
#               ../extracted_data/{sim}_hedgehog.pkl

python plot_ic_figures.py                # fig_ics.pdf, 
                                         # fig_amr_onecol.pdf
python plot_sim_showoff3_alt.py          # fig_face-comp_slices_alt.pdf
python plot_disk_mass_evolution.py       # fig_mass-ev.pdf, 
                                         # fig_mass-ev_late.pdf,
                                         # fig_net-gas-mass.pdf
python plot_sim_showoff1.py              # fig_edge-ev_fid.pdf
python plot_sim_showoff2.py              # fig_edge-comp.pdf
python plot_tcool_dist.py                # fig_tcool-mass-dist_cumm-big.pdf,
                                         # fig_tcool-mass-dis_tctfff-var.pdf
python plot_hedgehog.py                  # fig_entropy-pressure-tctff-mass_fid.pdf,
                                         # fig_CGM90_entropy-tctff-vel.pdf

