set(GEOXY_FLOCQ_RELEASES https://gitlab.inria.fr/flocq/flocq/-/archive)

FetchContent_Declare(
  flocq
  URL ${GEOXY_FLOCQ_RELEASES}/flocq-4.2.2/flocq-flocq-4.2.2.tar.gz
  URL_HASH
    SHA256=3ef17117828f4596c8c7e0dab4f5f3f1d0afa01e00b840ccbedff95632092040
)

FetchContent_MakeAvailable(flocq)
