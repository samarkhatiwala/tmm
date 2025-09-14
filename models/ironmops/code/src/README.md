iMOPS - MOPS with a simple iron cycle 
This in an unpublished version, which may still undergo some changes.

Iris Kriest, 20 August 2025, ikriest@geomar.de

Main changes to MOPS: 

- 2 additional tracers (dissolved and particulate iron: dfe and pfe), requires files for initialisation (dfeini.petsc and pfeini.petsc)
- supply of dFe from atmospheric and hydrothermal sources (requires files Fes_*)
- biogeochemical parameters and cycling (see BGC_INI.F and BGC_MODEL.F)
- additional cost function component (dFe) for optimisation (requires files dfe_obs.petsc and dfe_weight.petsc) 

Short desciption:

MOPS was extended by a simple iron component that simulates the cycling of dissolved and particulate iron (dFe and pFe, respectively). 
The iron components are calculated in units of μmol Fe m−3. The interaction of iron with plankton components is parameterised by a simplified version of the model by Nickelsen et al. (2015, https://doi.org/10.5194/gmd-8-1357-2015); here, we skipped their iron-dependent C:Chl ratio and modified the iron-dependent light affinity of phytoplankton as described in detail elsewhere. 
We assume that dissolved organic matter does not contain any iron - hence, the release or consumption of this component does not affect iron cycling. For the plankton components we assume a constant iron-to-phosphorus ratio. 
This constant proportion propagates into the detritus components; however, scavenging of dissolved iron can increase the Fe:P ratio in detritus. 
Similar to earlier model studies we assume that dissolved iron can be released from its particulate counterpart (Aumont et al., 2015, https://doi.org/10.5194/gmd-8-2465-2015; Nickelsen et al., 2015; Somes et al., 2021, https://doi.org/10.1029/2021GB006948), but not beyond the ratio assigned to plankton.
Scavenging and precipitation of dissolved iron are parameterised as in Tagliabue and Völker (2011, https://doi.org/10.5194/bg-8-3025-2011) and Aumont et al. (2015), including their approach of a variable ligand concentration. 
Iron is supplied at the model’s boundaries (surface and seafloor) via two external sources, namely atmospheric deposition and hydrothermal release, following Somes et al. (2021). 
It is exported from the model domain through burial of particulate iron at the sea floor. 

The file BGC_INI.F contains some first-guess default parameters.
The model's parameters (the Fe:P ratio, DOP-dependence of ligands, scavenging parameters) were also optimised against observations of inorganic nutrients, oxygen and dissolved iron; if you want to use those contact ikriest@geomar.de.
