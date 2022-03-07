#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.chunks = 1

IONICE = 'ionice -c2 -n7'
