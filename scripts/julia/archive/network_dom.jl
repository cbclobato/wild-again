using FlashWeave
#data_path = "outputs/r/311_dif_ab/ab_domc50_rep_all.csv"
#meta_data_path = "outputs/r/311_dif_ab/metadata_domc50_rep_all.csv"
#netw_results = learn_network(
#    data_path,
#    meta_data_path,
#    normalize=true,
#    sensitive=true,
#    alpha=0.99,
#    heterogeneous=false,
#    track_rejections=false,
#    n_obs_min=0,
#    verbose=true
#    )
#output_path = "outputs/flashweave/network_domc50_rep_all.gml"
#save_network(output_path, netw_results)

data_path = "C:\\\\Users\\pongo\\Documents\\GitHub\\cannabis-seed-microbiome\\outputs\\r\\311_dif_ab\\ab_domc50_rep_all.csv"
meta_data_path = "C:\\\\Users\\pongo\\Documents\\GitHub\\cannabis-seed-microbiome\\outputs\\r\\311_dif_ab\\metadata_domc50_rep_all.csv"
netw_results = learn_network(
    data_path,
    meta_data_path,
    normalize=true,
    sensitive=true,
    alpha=0.99,
    heterogeneous=false,
    track_rejections=false,
    n_obs_min=0,
    verbose=true
    ) #   prec=64
output_path = "C:\\\\Users\\pongo\\Documents\\GitHub\\cannabis-seed-microbiome\\outputs\\flashweave\\network_domc50_rep_all.gml"
save_network(output_path, netw_results)