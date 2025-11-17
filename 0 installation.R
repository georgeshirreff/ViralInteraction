# Epidemiological model of interactions between two seasonal respiratory viruses
# Copyright (C) 2025 George Shirreff
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.



# installation

# install packages
install.packages("tidyverse") 
install.packages("odin")
install.packages("FME")

install.packages("stringi")
install.packages("ggh4x")
install.packages("coda")
install.packages("data.table")
install.packages("ggstance")
install.packages("ggpubr")

# install the odin model
# to check that this works

# can take a few minutes
source("odin_m5.R")