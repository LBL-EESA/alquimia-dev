include(add_alquimia_batch_chem_benchmark)

# NOTE: for the moment, the batch chemistry driver doesn't support geochemical 
# NOTE: conditions specified in .cfg files, so we have disabled tests of this
# NOTE: capability (*_ac_*.cfg).

if (ALQUIMIA_HAVE_PFLOTRAN)
#  add_alquimia_batch_chem_benchmark(calcite_short_ac_pflotran calcite-short-ac-pflotran.cfg)
  add_alquimia_batch_chem_benchmark(calcite_short_pc_pflotran calcite-short-pc-pflotran.cfg)
  add_alquimia_batch_chem_benchmark(calcite_vf_pc_pflotran calcite-vf-pc-pflotran.cfg)
#  add_alquimia_batch_chem_benchmark(general_reaction_ac_pflotran general-reaction-ac-pflotran.cfg)
  add_alquimia_batch_chem_benchmark(general_reaction_pc_pflotran general-reaction-pc-pflotran.cfg)
#  add_alquimia_batch_chem_benchmark(ion_exchange_valocchi_ac_pflotran ion-exchange-valocchi-ac-pflotran.cfg)
  add_alquimia_batch_chem_benchmark(ion_exchange_valocchi_pc_pflotran ion-exchange-valocchi-pc-pflotran.cfg)
#  add_alquimia_batch_chem_benchmark(isotherms_ac_pflotran isotherms-ac-pflotran.cfg)
  add_alquimia_batch_chem_benchmark(isotherms_pc_pflotran isotherms-pc-pflotran.cfg)
#  add_alquimia_batch_chem_benchmark(surface_complexation_2_ac_pflotran surface-complexation-2-ac-pflotran.cfg)
  add_alquimia_batch_chem_benchmark(surface_complexation_2_pc_pflotran surface-complexation-2-pc-pflotran.cfg)
  add_alquimia_batch_chem_benchmark(microbial_pflotran microbial-pflotran.cfg)
endif()

if (ALQUIMIA_HAVE_CRUNCHFLOW)
#  add_alquimia_batch_chem_benchmark(calcite_short_ac_crunch calcite-short-ac-crunch.cfg)
  add_alquimia_batch_chem_benchmark(calcite_short_cc_crunch calcite-short-cc-crunch.cfg)
#  add_alquimia_batch_chem_benchmark(calcite_vf_ac_crunch calcite-vf-ac-crunch.cfg)
  add_alquimia_batch_chem_benchmark(calcite_vf_cc_crunch calcite-vf-cc-crunch.cfg)
#  add_alquimia_batch_chem_benchmark(general_reaction_ac_crunch general-reaction-ac-crunch.cfg)
  add_alquimia_batch_chem_benchmark(general_reaction_cc_crunch general-reaction-cc-crunch.cfg)
#  add_alquimia_batch_chem_benchmark(ion_exchange_valocchi_ac_crunch ion-exchange-valocchi-ac-crunch.cfg)
  add_alquimia_batch_chem_benchmark(ion_exchange_valocchi_cc_crunch ion-exchange-valocchi-cc-crunch.cfg)
#  add_alquimia_batch_chem_benchmark(isotherms_ac_crunch isotherms-ac-crunch.cfg)
  add_alquimia_batch_chem_benchmark(isotherms_cc_crunch isotherms-cc-crunch.cfg)
#  add_alquimia_batch_chem_benchmark(surface_complexation_2_ac_crunch surface-complexation-2-ac-crunch.cfg)
  add_alquimia_batch_chem_benchmark(surface_complexation_2_cc_crunch surface-complexation-2-cc-crunch.cfg)
endif()

