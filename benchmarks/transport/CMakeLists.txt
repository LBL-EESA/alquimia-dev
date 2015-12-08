# NOTE: Currently all the benchmarks are just smoke tests.

include(add_alquimia_transport_benchmark)

if (ALQUIMIA_HAVE_PFLOTRAN)
  add_alquimia_transport_benchmark(tracer_1d_pflotran tracer-1d-pflotran.cfg)
endif()

if (ALQUIMIA_HAVE_CRUNCHFLOW)
endif()

