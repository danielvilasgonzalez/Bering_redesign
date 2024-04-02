# Adapting monitoring to a changing seascape: increasing the efficiency, flexibility, and continuity of bottom trawl surveys in the Bering Sea and beyond

## Description
This repository provides the code used for an In Prep manuscript by Daniel Vilas, Lewis A.K. Barnett, André E. Punt, Zack S. Oyafuso, Lukas DeFilippo, Margaret C. Siple, Leah S. Zacher, and Stan Kotwicki entitled "**Which traditional sampling design best estimates abundance across species in a rapidly changing ecosystem?**".

## Material and Methods

We investigated whether defining survey boundaries based on historical and future environmental conditions improves the precision and accuracy of abundance estimates in a multispecies survey. We fitted univariate spatiotemporal species distribution models to 16 stocks (14 species) using historical observations of fishery-independent bottom trawl survey catch-per-unit-effort and sea bottom temperature in the eastern and northern Bering Sea from 1982 to 2022. We used spatiotemporal models to simulate historical and future survey data from these models and optimize stratuma boundaries and sample allocation for abundance estimation under a variety of environmental conditions. We then compared simulated abundance estimates to the simulated true abundance among sampling designs and future temperature scenarios.

## Species/Stocks Included

The species and stocks set included in the manuscript are 10 groundfish and 4 crab species (6 stocks) of the Bering Sea:

| Stock                               | Scientific Name                     | Common Name                           |
|-------------------------------------|-------------------------------------|---------------------------------------|
| arrowtooth flounder                 | *Atheresthes stomias*               | arrowtooth flounder                   |
| Arctic cod                          | *Boreogadus saida*                  | Arctic cod                            |
| Tanner crab                         | *Chionoecetes bairdi*               | Tanner crab                           |
| snow crab                           | *Chionoecetes opilio*               | snow crab                             |
| saffron cod                         | *Eleginus gracilis*                 | saffron cod                           |
| Alaska pollock                      | *Gadus chalcogrammus*               | Alaska pollock                        |
| Pacific cod                         | *Gadus macrocephalus*               | Pacific cod                           |
| flathead sole                       | *Hippoglossoides elassodon*         | flathead sole                         |
| Bering flounder                     | *Hippoglossoides robustus*          | Bering flounder                       |
| northern rock sole                  | *Lepidopsetta polyxystra*           | northern rock sole                    |
| yellowfin sole                      | *Limanda aspera*                    | yellowfin sole                        |
| Pribilof Islands blue king crab     | *Paralithodes platypus*             | blue king crab                        |
| St. Matthew Island blue king crab   | *Paralithodes platypus*             | blue king crab                        |
| Pribilof Islands red king crab      | *Paralithodes camtschaticus*        | red king crab                         |
| Bristol Bay red king crab           | *Paralithodes camtschaticus*        | red king crab                         |
| Alaska plaice                       | *Pleuronectes quadrituberculatus*   | Alaska plaice                         |


## Sampling designs

Stratification scheme and station allocation information for each sampling design. The “optimized” stratification schemes represent the multispecies optimal design. All sampling designs consist of 15 strata and 520 samples.


| Stratification scheme | Stratification factors                       | Sampling allocation       |
|-----------------------|----------------------------------------------|---------------------------|
| existing              | depth and geographical subregion             | fixed                     |
| existing              | depth and geographical subregion             | fixed                     |
| existing              | depth and geographical subregion             | fixed                     |
| optimized             | depth                                        | fixed                     |
| optimized             | depth                                        | fixed                     |
| optimized             | depth                                        | fixed                     |
| optimized             | variance of sea bottom temperature           | fixed                     |
| optimized             | variance of sea bottom temperature           | fixed                     |
| optimized             | variance of sea bottom temperature           | fixed                     |
| optimized             | depth and variance of sea bottom temperature | fixed                     |
| optimized             | depth and variance of sea bottom temperature | fixed                     |
| optimized             | depth and variance of sea bottom temperature | fixed                     |

