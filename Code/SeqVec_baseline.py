import numpy as np
from allennlp.commands.elmo import ElmoEmbedder
from pathlib import Path
import datetime
from multiprocessing import Pool

def get_embedding(sequence):
    return embedder.embed_sentence(sequence)

model_dir = Path('Path of pretrain SeqVec uniref50_v2 file')
weights = model_dir / 'weights.hdf5'
options = model_dir / 'options.json'

seq = np.load("sequence dataset .npy file")

embedder = ElmoEmbedder(options, weights)

seqvec_spike_embed_7k = []

start = datetime.datetime.now()
with Pool(16) as p:
    seqvec_spike_embed_7k = p.map(get_embedding, seq)
end = datetime.datetime.now()
print((end-start).total_seconds())

print('length of embeddings', len(seqvec_spike_embed_7k))
np.save(".npy file name to save embeddings", seqvec_spike_embed_7k)

print('done')

