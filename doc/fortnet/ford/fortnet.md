macro:
        DEBUG=0
        WITH_MPI=1
	WITH_SOCKETS=0
        RELEASE=0.1
preprocess: true
src_dir:
        ../../../prog/fortnet
output_dir: ./doc
project_github: https://github.com/vanderhe/fortnet
project_website: https://github.com/vanderhe/fortnet
summary: Fortnet: A Behler-Parrinello-Neural-Network implementation.
author: T. W. van der Heide
preprocessor: ../../../external/fypp/bin/fypp
include: ../../../prog/fortnet/include
predocmark: >
display: public
         protected
proc_internals:
        false
source: true
graph: true
search: false
license: by
warn: true
