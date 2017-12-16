# scpsdp_alm
Simple ALM solver for Burer-Monteiro approach on sparsest cut SDP

The basic flow is to use `run_sparsest_cut.m` to run one of the
algorithms implemented. See the first few lines of this script
for parameters you can/should change; in particular `graph_name`
can be set to any one of the graph file names in the
`data/graphs` directory, and `semimetric_mode` can be either
`lp`, `sdp`, or `bm`.

The other `.m` files included are mostly for utility puposes.
