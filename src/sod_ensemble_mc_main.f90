!*******************************************************************************
!    Copyright (c) 2025, Salvador R.G. Balestra
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!******************************************************************************

program sod_ensemble_mc
    use sod_ensemble_mc_mod, only: run_sod_ensemble_mc
    implicit none
    call run_sod_ensemble_mc()
end program sod_ensemble_mc
