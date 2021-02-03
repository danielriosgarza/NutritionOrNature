# NutritionOrNature

Repository of Figures and code used in the manuscript: [Nutrition or nature: disentangling the complex forces shaping prokaryote pan-genomes](https://www.biorxiv.org/content/10.1101/2020.12.14.422685v3)

Below is a guide to enable readers to reproduce our results:

1) The environment ball:

We simulated pan-reactomes in a random ball containing the relative concentration of metabolites. The environments we used in the manuscript are listed in [TableS3](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Table_S3.xlsx). We used the script: [generate_env_ball.py](xxx.xxx) 

2) Create family-level pan-reactomes:

We begin with a database of draft reconstractions of strain-level genome-scale metabolic models (GSMMs) listed in [TableS2](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Table_S2.xlsx). These contain strains that belong to families with 25 or more sequenced genomes that belog to different species. The models were reconstructed using [ModelSEED](https://modelseed.org/) and are available on request.

To merge models belonging to the same families, we used the script [create_ensemble_model.py](xxx.xxx). An example of a model reconstructed for the Aeromonadaceae family that is extensively described in the manuscript is available here: [Aeromonadaceae_reactome](xxxx.xxx)


3) The toy model example

4) Analyze the panEFMs

5) Get Elastic net predictions
