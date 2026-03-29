import numpy as np
import anndata as ad
from scipy import sparse
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


class H5ADStreamLoader:
    """Stream mini-batches from one or many .h5ad files without full in-memory concat.

    Notes
    -----
    - Uses ``anndata.read_h5ad(..., backed='r')`` so expression matrix stays on disk.
    - Emits NumPy arrays (or sparse CSR matrices) for model-side training loops.
    - Supports lightweight row filtering on obs columns.
    """

    def __init__(
        self,
        paths: Sequence[str],
        batch_size: int = 1024,
        obs_filter: Optional[Dict[str, Iterable[str]]] = None,
        dense_output: bool = True,
        shuffle_within_file: bool = True,
        seed: int = 0,
    ) -> None:
        if not paths:
            raise ValueError("`paths` cannot be empty.")
        if batch_size <= 0:
            raise ValueError("`batch_size` must be > 0.")

        self.paths = list(paths)
        self.batch_size = int(batch_size)
        self.obs_filter = obs_filter or {}
        self.dense_output = dense_output
        self.shuffle_within_file = shuffle_within_file
        self.rng = np.random.default_rng(seed)

    def _build_mask(self, obs) -> np.ndarray:
        mask = np.ones(obs.shape[0], dtype=bool)
        for key, allowed_values in self.obs_filter.items():
            if key not in obs.columns:
                raise KeyError(f"obs column `{key}` not found.")
            allowed_set = set(allowed_values)
            mask &= obs[key].isin(allowed_set).to_numpy()
        return mask

    def _slice_x(self, x):
        if self.dense_output:
            if sparse.issparse(x):
                return x.toarray()
            return np.asarray(x)
        if sparse.issparse(x):
            return x.tocsr()
        return sparse.csr_matrix(np.asarray(x))

    def batches(self, epochs: int = 1) -> Iterator[Tuple[np.ndarray, np.ndarray, str]]:
        """Yield (X_batch, obs_index_batch, source_path)."""
        if epochs <= 0:
            raise ValueError("`epochs` must be > 0.")

        for _ in range(epochs):
            file_order = np.arange(len(self.paths))
            self.rng.shuffle(file_order)

            for file_idx in file_order:
                path = self.paths[file_idx]
                adata = ad.read_h5ad(path, backed="r")

                try:
                    mask = self._build_mask(adata.obs)
                    indices = np.where(mask)[0]

                    if self.shuffle_within_file:
                        self.rng.shuffle(indices)

                    for start in range(0, len(indices), self.batch_size):
                        batch_idx = indices[start : start + self.batch_size]
                        if len(batch_idx) == 0:
                            continue
                        x_batch = adata.X[batch_idx]
                        yield self._slice_x(x_batch), batch_idx, path
                finally:
                    adata.file.close()

    def iter_obs_chunks(
        self,
        columns: Sequence[str],
        chunk_size: int = 50000,
    ) -> Iterator[Tuple[Dict[str, np.ndarray], str]]:
        """Stream selected obs columns in chunks.

        Useful for computing dataset-level statistics without loading X.
        """
        if chunk_size <= 0:
            raise ValueError("`chunk_size` must be > 0.")

        for path in self.paths:
            adata = ad.read_h5ad(path, backed="r")
            try:
                for col in columns:
                    if col not in adata.obs.columns:
                        raise KeyError(f"obs column `{col}` not found in {path}.")

                n = adata.n_obs
                for start in range(0, n, chunk_size):
                    end = min(start + chunk_size, n)
                    block = {
                        col: adata.obs[col].iloc[start:end].to_numpy()
                        for col in columns
                    }
                    yield block, path
            finally:
                adata.file.close()
