from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys



def plot_genome(path_to_gb):
  """
  Takes a path to the genbank file.
  Plots the genome in given gb file.

  This function is based on an example given in the biopython tutorial/cookbook.
  See http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc264 for original example.
  """
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
      gd_feature_set.add_feature(feature, sigil="BOX", 
                                 color="brown", line_width=30, strand=None)


  gd_diagram.draw(format="circular", circular=True, pagesize=(10*cm,10*cm),
                  start=0, end=len(record), circle_core = .8)
  gd_diagram.write(path_to_gb[:-2]+"pdf", "PDF")


if __name__ == '__main__':
  plot_genome(sys.argv[1])