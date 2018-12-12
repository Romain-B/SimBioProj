from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys



def plot_genome(path_to_gb):
  record = SeqIO.read(path_to_gb, "genbank")

  gd_diagram = GenomeDiagram.Diagram(record.id)
  gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
  gd_feature_set = gd_track_for_features.new_set()

  for feature in record.features:
    if feature.type == "gene":
          #Exclude this feature
          #continue
      if len(gd_feature_set) % 2 == 0:
          color = colors.blue
      else:
          color = colors.lightblue
      gd_feature_set.add_feature(feature, sigil="ARROW", arrowshaft_height=1.0, arrowshaft_width=.5,
                                 color=color, label=True,
                                 label_size = 14, label_angle=90 )


  #for feature in record.features:
    if feature.type == "barrier":
          #Exclude this feature
          #continue
      gd_feature_set.add_feature(feature, sigil="BOX", 
                                 color="brown", line_width=30, strand=None)


  gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm),
                  start=0, end=len(record), circle_core = .8)
  gd_diagram.write(path_to_gb[:-2]+"pdf", "PDF")
#gd_diagram.write("plasmid_circular_nice.eps", "EPS")
#gd_diagram.write("plasmid_circular_nice.svg", "SVG")

if __name__ == '__main__':
  plot_genome(sys.argv[1])