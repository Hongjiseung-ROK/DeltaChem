import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing, global_add_pool

class EquivariantLayer(MessagePassing):
    def __init__(self, in_features, out_features):
        super(EquivariantLayer, self).__init__(aggr='add')
        self.msg_mlp = nn.Sequential(
            nn.Linear(in_features * 2 + 1, out_features),
            nn.Silu(),
            nn.Linear(out_features, out_features)
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(in_features + out_features, out_features),
            nn.Silu(),
            nn.Linear(out_features, out_features)
        )

    def forward(self, x, edge_index, edge_attr):
        return self.propagate(edge_index, x=x, edge_attr=edge_attr)

    def message(self, x_i, x_j, edge_attr):
        # edge_attr is squared distance
        msg = torch.cat([x_i, x_j, edge_attr], dim=-1)
        return self.msg_mlp(msg)

    def update(self, aggr_out, x):
        new_x = torch.cat([x, aggr_out], dim=-1)
        return self.node_mlp(new_x)

class EGNN(nn.Module):
    def __init__(self, num_features=1, hidden_channels=64, num_layers=4, dropout=0.1):
        super(EGNN, self).__init__()
        self.embedding = nn.Linear(num_features, hidden_channels)
        self.layers = nn.ModuleList([
            EquivariantLayer(hidden_channels, hidden_channels) for _ in range(num_layers)
        ])
        self.dropout = nn.Dropout(p=dropout)
        self.output = nn.Sequential(
            nn.Linear(hidden_channels, 32),
            nn.Silu(),
            nn.Linear(32, 1) # Predicting Delta-Energy
        )

    def forward(self, data, training=True):
        x, edge_index, pos, batch = data.x, data.edge_index, data.pos, data.batch
        
        # Calculate edge_attr as squared distances
        row, col = edge_index
        dist_sq = torch.sum((pos[row] - pos[col])**2, dim=-1, keepdim=True)
        
        x = self.embedding(x)
        for layer in self.layers:
            x = layer(x, edge_index, dist_sq)
            if training:
                x = self.dropout(x)
        
        x = global_add_pool(x, batch)
        return self.output(x)

    def predict_with_uncertainty(self, data, n_samples=20):
        """MC-Dropout for uncertainty estimation."""
        self.train() # Enable dropout
        preds = [self.forward(data, training=True) for _ in range(n_samples)]
        preds = torch.stack(preds)
        mean = preds.mean(dim=0)
        std = preds.std(dim=0)
        return mean, std
