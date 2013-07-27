#!/usr/bin/env python
import os
import lagrange
data = """\
### begin data
{'area_adjacency': [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]],
 'area_dispersal': [[[1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0]],
                    [[1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0]],
                    [[1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0]],
                    [[1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0]],
                    [[1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0]]],
 'area_labels': ['K', 'O', 'M', 'H'],
 'base_rates': '__estimate__',
 'dispersal_durations': [0.5,
                         1.3999999999999999,
                         1.8000000000000003,
                         1.3999999999999995,
                         4.9000000000000004],
 'dm_symmetric_entry': True,
 'excluded_ranges': [],
 'lagrange_version': '20130526',
 'max_range_size': 2,
 'model_name': 'psychotria_M1_strat.lagrange',
 'newick_trees': [{'included': '__all__',
                   'name': 'Tree0',
                   'newick': '((((((((P_hawaiiensis_WaikamoiL1:0.9656850499,P_mauiensis_Eke:0.9656850499):0.7086257935,(P_fauriei2:1.230218511,P_hathewayi_1:1.230218511):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.185759997):0.302349378,(P_greenwelliae07:1.131363255,P_greenwelliae907:1.131363255):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.994011054,P_hawaiiensis_Makaopuhi:1.994011054):0.7328279804,P_mariniana_Oahu:2.726839034):0.2574151709,P_mariniana_Kokee2:2.984254205):0.4601084855,P_wawraeDL7428:3.444362691):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.479300491,P_hobdyi_Kuia:2.479300491):2.432497733):0.2873119899,((P_hexandra_K1:2.363984189,P_hexandra_M:2.363984189):0.4630447802,P_hexandra_Oahu:2.826939991):2.372081244);',
                   'root_age': 5.2}],
 'ranges': [(),
            (0,),
            (0, 1),
            (0, 2),
            (0, 3),
            (1,),
            (1, 2),
            (1, 3),
            (2,),
            (2, 3),
            (3,)],
 'taxa': ['P_mariniana_Kokee2',
          'P_mariniana_Oahu',
          'P_mariniana_MauiNui',
          'P_hawaiiensis_Makaopuhi',
          'P_wawraeDL7428',
          'P_kaduana_PuuKukuiAS',
          'P_mauiensis_PepeAS',
          'P_hawaiiensis_WaikamoiL1',
          'P_mauiensis_Eke',
          'P_fauriei2',
          'P_hathewayi_1',
          'P_kaduana_HawaiiLoa',
          'P_greenwelliae07',
          'P_greenwelliae907',
          'P_grandiflora_Kal2',
          'P_hobdyi_Kuia',
          'P_hexandra_K1',
          'P_hexandra_M',
          'P_hexandra_Oahu'],
 'taxon_range_data': {'P_fauriei2': (1,),
                      'P_grandiflora_Kal2': (0,),
                      'P_greenwelliae07': (0,),
                      'P_greenwelliae907': (0,),
                      'P_hathewayi_1': (1,),
                      'P_hawaiiensis_Makaopuhi': (3,),
                      'P_hawaiiensis_WaikamoiL1': (2,),
                      'P_hexandra_K1': (0,),
                      'P_hexandra_M': (0,),
                      'P_hexandra_Oahu': (1,),
                      'P_hobdyi_Kuia': (0,),
                      'P_kaduana_HawaiiLoa': (1,),
                      'P_kaduana_PuuKukuiAS': (2,),
                      'P_mariniana_Kokee2': (0,),
                      'P_mariniana_MauiNui': (2,),
                      'P_mariniana_Oahu': (1,),
                      'P_mauiensis_Eke': (2,),
                      'P_mauiensis_PepeAS': (2,),
                      'P_wawraeDL7428': (0,)}}
### end data
"""

i = 0
while 1:
    if not i:
        outfname = "psychotria_M1_strat.results.txt"
    else:
        outfname = "psychotria_M1_strat.results-"+str(i)+".txt"
    if not os.path.exists(outfname): break
    i += 1
outfile = open(outfname, "w")
lagrange.output.log(lagrange.msg, outfile, tee=True)
model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(data)
lagrange.output.ascii_tree(outfile, tree, model, data, tee=True)
if base_rates != "__estimate__":
    d, e = base_rates
else:
    d, e = lagrange.output.optimize_dispersal_extinction(outfile, tree, model, tee=True)
if nodelabels:
    if nodelabels == "__all__":
        nodelabels = None
    lagrange.output.ancsplits(outfile, tree, model, d, e, nodelabels=nodelabels, tee=True)
