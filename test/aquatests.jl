import Aqua

using RenewalDiD

Aqua.test_all(RenewalDiD)
Aqua.test_undocumented_names(RenewalDiD; broken=true)
