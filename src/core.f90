!*******************************************************************************
! Shared SOD core API: legacy names on top of modern reusable modules.
!*******************************************************************************
! Module `core` reexports the shared high-level API used by the modern SOD workflows.
! El módulo `core` reexporta la API compartida de alto nivel usada por los flujos modernos de SOD.
module core
    use consts, only: dp, ip, kB_eVk
    use utils, only: binomial_int64, next_combination, random_subset, format_level_directory, &
        join_path, ensure_directory_exists
    use inputs, only: sgo_file_data, insod_file_data, read_sgo_file, read_insod_file
    use configurations, only: canonicalize_subset, collect_subset_stabilizer_operators, &
        enumerate_unique_subsets, find_subset_index, read_outsod_file, read_outsod_header, write_outsod_file, write_outsod_unit
    use structure_io, only: write_vasp_configuration
    use symmetry, only: symmetry_initialize, symmetry_finalize, symmetry_get_matrix, &
        symmetry_generate_matrix_from_files, symmetry_restrict_supercell_operators, symmetry_validate_matrix, &
        symmetry_wrap_fractional_vector
    implicit none
    public
end module core
