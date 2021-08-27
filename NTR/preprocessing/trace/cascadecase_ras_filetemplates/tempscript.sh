#!/bin/bash
cd /bigwork/nhkcmany/12_gwk_tracecase
/sw/TRACE/openmpi-1.3.1/bin/mpiexec /sw/TRACE/trace_7.1.22_openmpi_MTU/TRACE -cgns /bigwork/nhkcmany/input/TRACE.cgns -lb /bigwork/nhkcmany/input/BALANCE_1PROC -o TRACE.lst. -pc
