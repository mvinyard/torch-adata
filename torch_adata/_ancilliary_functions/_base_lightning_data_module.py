
from pytorch_lightning import LightningDataModule
from torch.utils.data import DataLoader
import os


class BaseLightningDataModule(LightningDataModule):
    def __init__(self, batch_size=2000, num_workers=os.cpu_count()):
        super().__init__()

        self.batch_size = batch_size
        self.num_workers = num_workers

    def train_dataloader(self):
        return DataLoader(
            self.train_dataset, batch_size=self.batch_size, num_workers=self.num_workers
        )

    def val_dataloader(self):
        return DataLoader(
            self.val_dataset, batch_size=self.batch_size, num_workers=self.num_workers
        )

    def test_dataloader(self):
        return DataLoader(
            self.test_dataset, batch_size=self.batch_size, num_workers=self.num_workers
        )

    def predict_dataloader(self):
        return DataLoader(
            self.dataset, batch_size=self.batch_size, num_workers=self.num_workers
        )