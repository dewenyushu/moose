#!/bin/bash

mpiexec -np 1 ./../../../combined-opt -i iaea_vp4_xyz.i Problem/restore_original_nonzero_pattern=false Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-no-reset-1procs.log
mpiexec -np 1 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-reset-1procs.log
mpiexec -np 1 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=true Problem/error_on_jacobian_nonzero_reallocation=true   > hash-1procs.log

mpiexec -np 2 ./../../../combined-opt -i iaea_vp4_xyz.i Problem/restore_original_nonzero_pattern=false Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-no-reset-2procs.log
mpiexec -np 2 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-reset-2procs.log
mpiexec -np 2 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=true Problem/error_on_jacobian_nonzero_reallocation=true   > hash-2procs.log

mpiexec -np 4 ./../../../combined-opt -i iaea_vp4_xyz.i Problem/restore_original_nonzero_pattern=false Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-no-reset-4procs.log
mpiexec -np 4 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-reset-4procs.log
mpiexec -np 4 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=true Problem/error_on_jacobian_nonzero_reallocation=true   > hash-4procs.log

mpiexec -np 8 ./../../../combined-opt -i iaea_vp4_xyz.i Problem/restore_original_nonzero_pattern=false Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-no-reset-8procs.log
mpiexec -np 8 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-reset-8procs.log
mpiexec -np 8 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=true Problem/error_on_jacobian_nonzero_reallocation=true   > hash-8procs.log

mpiexec -np 16 ./../../../combined-opt -i iaea_vp4_xyz.i Problem/restore_original_nonzero_pattern=false Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-no-reset-16procs.log
mpiexec -np 16 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=false Problem/error_on_jacobian_nonzero_reallocation=false   > preallocation-reset-16procs.log
mpiexec -np 16 ./../../../combined-opt -i iaea_vp4_xyz.i  Problem/restore_original_nonzero_pattern=true Problem/use_hash_table_matrix_assembly=true Problem/error_on_jacobian_nonzero_reallocation=true   > hash-16procs.log
