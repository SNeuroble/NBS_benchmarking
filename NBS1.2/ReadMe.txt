1. Unzip NBS1.x.zip
2. Launch Matlab 
3. Add folder NBS1.x to your Matlab path (File->Set Path->Add Folder)
4. Type NBS at the Matlab prompt

Changes to NBS v1.2
- Fixed the NaN bug in NBSglm that arises when the dependent variable is identically zero for all observations. Force any test statistic value that is NaN to zero. This occurs when the connectivity matrix contains a cell that is identically zero for all observations. Future versions will allow a preliminary test to be performed to exclude from statistical testing any elements in the connectivity matrix that are identically (or sufficiently close to) zero for all observations.  
- Fixed the p-value bug in NBSfdr. Previously, the p-value counter was not incremented if the observed test statistic was equal to a randomized test statistic. 
- Revised NBSstats to enable fractional component sizes to be shown, which is possible if component size is measured with the intensity option. 

New features to NBS v1.1:
- Graphical user interface
- Statistical design specified with general linear model. Covers most statistical designs, including t-test, F-test, ANOVA, ANCOVA, multiple linear regression, etc., and enables modeling of nuisance covariates. 
- Exchange blocks for constraining permutations in repeated measures designs. 
- NBSview, a basic network viewer modeled on SPMresults. 
- Option to measure netowrk size with intensity.
- False discovery rate (FDR) based on nonparametric p-values.
- Permutations stored enabling fast re-computation when experimenting with different test-statistic thresholds. 

 Copyright (C) 2012

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

azalesky@unimelb.edu.au 
