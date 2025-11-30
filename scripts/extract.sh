#!/bin/bash
#*******************************************************************************
#    Copyright (c) 2025, Salvador R.G. Balestra
#
#    This file is part of the SOD package.
#
#    SOD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#******************************************************************************

# Collects the final lattice energy from every *.vasp.gout file and writes them to ENERGIES.
for f in *.vasp.gout; do
	# skip files that indicate non-converged minimum
	if grep -qF 'Conditions for a minimum have not been satisfied.' "$f"; then
	    echo "N/A # $f"
		continue
	fi
	x=$(grep 'l e' "$f" | tail -n1 | awk '{print $4}')
	echo "$x # $f"
done > ENERGIES
