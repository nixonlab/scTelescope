#! /usr/bin/env python

import os


bamfile = "Esta.Chingadera.Tiene.Unchingo.depuntos.bam"
data_name = '.'.join(os.path.basename(bamfile).split('.')[0:-1])
extension = os.path.basename(bamfile).split('.')[-1]
print(data_name)
print(extension)
exit()
