#!/usr/bin/env nextflow

/*
 * process
 */
process genericProcess {

    container 'singlecellopenproblems/openproblems'
    queue 'default-2a6498e0-72f2-11eb-b638-06d5c9fd5439'

    """
    echo hello world
    """
}
