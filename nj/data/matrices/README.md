# Empirical protein rate matrices

Exchangeabilities and equilibrium frequencies for the built-in protein models, in PAML layout: the strict
lower triangle of the symmetric exchangeability matrix by rows (19 rows of 1, 2, …, 19 numbers), then the 20
amino-acid frequencies, in the order `A R N D C Q E G H I L K M F P S T W Y V`. The files are unmodified
copies of the ones listed on the page below; `ReversibleModel::from_paml` scales them to one substitution per
unit time, so only the ratios of the numbers matter.

Retrieved 2026-11-26 from
<http://giphy.pasteur.fr/empirical-models-of-amino-acid-substitution/> (page "last update: 22.09.30"), which
lists, for each model, the reference and where the data came from. The `Data source` column repeats what that
page states; nj.rs did not fetch from those sources itself.

| nj.rs model (CLI name) | File | Reference | Data source, as stated on the page |
|---|---|---|---|
| `Lg` (`lg`) | `LG.dat` | Le SQ, Gascuel O (2008) An improved general amino acid replacement matrix. *Mol Biol Evol* 25:1307–1320. doi:[10.1093/molbev/msn067](https://doi.org/10.1093/molbev/msn067) | <http://www.atgc-montpellier.fr/download/datasets/models/lg_LG.PAML.txt> |
| `Wag` (`wag`) | `WAG.dat` | Whelan S, Goldman N (2001) A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. *Mol Biol Evol* 18:691–699. doi:[10.1093/oxfordjournals.molbev.a003851](https://doi.org/10.1093/oxfordjournals.molbev.a003851) | <https://www.ebi.ac.uk/goldman-srv/WAG/wag.dat> |
| `Jtt` (`jtt`) | `JTT.dat` | Jones DT, Taylor WR, Thornton JM (1992) The rapid generation of mutation data matrices from protein sequences. *Comput Appl Biosci* 8:275–282. doi:[10.1093/bioinformatics/8.3.275](https://doi.org/10.1093/bioinformatics/8.3.275) | MrBayes `src/model.c` (`/* jones */`): <https://raw.githubusercontent.com/NBISweden/MrBayes/master/src/model.c> |
| `JttDcmut` (`jtt-dcmut`) | `DCMut-JTT.dat` | Kosiol C, Goldman N (2005) Different versions of the Dayhoff rate matrix. *Mol Biol Evol* 22:193–199. doi:[10.1093/molbev/msi005](https://doi.org/10.1093/molbev/msi005) | <https://www.ebi.ac.uk/goldman-srv/dayhoff/jtt-dcmut.dat> |
| `Dayhoff` (`dayhoff`) | `Dayhoff.dat` | Dayhoff MO, Schwartz RM, Orcutt BC (1978) A model of evolutionary change in proteins. In: Dayhoff MO (ed.), *Atlas of Protein Sequence and Structure*, National Biomedical Research Foundation, Washington DC, Vol. 5, pp. 345–352 | <https://www.ebi.ac.uk/goldman-srv/dayhoff/dayhoff-paml.dat> |
| `DayhoffDcmut` (`dayhoff-dcmut`) | `DCMut-Dayhoff.dat` | Kosiol C, Goldman N (2005), as above | <https://www.ebi.ac.uk/goldman-srv/dayhoff/dayhoff-dcmut.dat> |

SHA-256 of the files as embedded:

```
8ffb9c24d30e5b332df96a53fd421258aa89092a2b9876c01ce4a9b0fb56d1fe  LG.dat
eed02827db884ef9bf52bf4a6f73df8ab8f0d5f1fe50fa8374f9c2bcd8cc2bde  WAG.dat
30ea3e3423d0bc23135a23d692a34a940489972b93010e5dd0ceebb918e3589f  JTT.dat
c4ddb19f52286f20c0b3fb054cee41500e56a511c8aa57e7d4c37a500135bc73  DCMut-JTT.dat
37a0b5e17fef50c1340f8d4b0369d3898e7a183e1ebeb6607088b005896a3dc3  Dayhoff.dat
2229cc9b9c52cca612a1da167204d590bee3d4714fb9fd3f86e365be5b0de024  DCMut-Dayhoff.dat
```

## Which version to use

Kosiol and Goldman (2005) found errors and rounding in the widely circulated Dayhoff and JTT rate matrices and
published recomputed versions ("DCMut"), which they recommend; `Dayhoff` and `Jtt` are kept because published
analyses used them. PAML's own `jones.dat` is not identical to the 1992 paper: it is based on an updated data set
from D. T. Jones, and the JTT file here comes from MrBayes rather than from the authors.

## Terms

The source page states no licence or terms of use for the files, and neither do the files. They are numerical
tables from the publications above; please cite those papers when you use the models. Whether redistribution inside
an MIT-licensed crate needs more than that, in particular for `JTT.dat`, whose stated source is the GPL-licensed
MrBayes, is for the maintainer to decide; this note is not legal advice.
