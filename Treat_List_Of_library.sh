#!/bin/bash
parallel --jobs 32 ./27.pe.pipeline.sh {1} :::: $1
