#!/usr/bin/env nextflow

/*
========================================================================================
    PHINDER: PHage Isolate characterizatioN, Discovery & Evaluation Resource
========================================================================================
    Github : https://github.com/tdoerks/PHINDER
    Author : Tyler Doerksen
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { PHINDER_PIPELINE } from './workflows/phinder_pipeline'

//
// WORKFLOW: Run main PHINDER analysis pipeline
//
workflow {
    PHINDER_PIPELINE()
}

/*
========================================================================================
    THE END
========================================================================================
*/
