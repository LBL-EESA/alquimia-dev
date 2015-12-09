include(add_alquimia_transport_benchmark)

if (ALQUIMIA_HAVE_PFLOTRAN)
  add_alquimia_transport_benchmark(tracer_1d_pflotran tracer-1d-pflotran.cfg)
  add_alquimia_transport_benchmark(tritium_1d_pflotran tritium-1d-pflotran.cfg)
  add_alquimia_transport_benchmark(calcite_1d_pflotran calcite-1d-pflotran.cfg)
  add_alquimia_transport_benchmark(isotherms_1d_pflotran isotherms-1d-pflotran.cfg)
  # The ion exchange below hangs during the reaction solve(!).
  #add_alquimia_transport_benchmark(ion_exchange_1d_pflotran ion-exchange-1d-pflotran.cfg)
  add_alquimia_transport_benchmark(farea_1d_pflotran farea-1d-pflotran.cfg)
endif()

if (ALQUIMIA_HAVE_CRUNCHFLOW)
endif()

