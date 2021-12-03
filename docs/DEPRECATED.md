### Deprecated functions
#### deprecated functions in module gt_file
1. `check_file_empty`, use `check_file_exist` instead.
2. `get_seqfile_prefix`, use `get_file_prefix` instead.
3. `perfect_open`, use `open_file`
4. `close_file`
5. `find_files`, use `find_file`

#### deprecated functions in module gt_path
1. `check_path_empty`, use `check_path_exist`
2. `abs_path`
3. `find_db_path`, use `find_db`

#### deprecated functions in module io_seq
1. `gc_to_dict`, use `seq2dict`
2. `seq_to_dict`, use `seq2dict`
3. `evaluate_genome`, use `assess_genome`

#### deprecated functions in module alter_seq
1. `select_seq`, use `select_seq_len`
2. `reorder_PE_fq`, use `sort_pe_fq`