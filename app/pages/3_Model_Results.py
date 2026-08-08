"""Model Results page: AUROC/AUPRC curves, confusion matrix, comparison
table.

Expects `data/processed/model_metrics.json`, written by the notebook's
Section 4 (Evaluate the Model), shaped as:

    {
      "models": {
        "<model name>": {
          "accuracy": float, "precision": float, "recall": float,
          "f1": float,
          "auroc": float, "fpr": [...], "tpr": [...],
          "auprc": float, "precision_curve": [...], "recall_curve": [...],
          "confusion_matrix": [[tn, fp], [fn, tp]]
        },
        ...
      }
    }
"""

import json
import sys
from pathlib import Path

import streamlit as st
import plotly.graph_objects as go

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

METRICS_PATH = REPO_ROOT / "data" / "processed" / "model_metrics.json"

st.set_page_config(page_title="PPI Prediction | Model Results", page_icon="🧬")
st.title("Model Results")

if not METRICS_PATH.exists():
    st.info(
        "Not generated yet. Once the notebook's Section 4 (Evaluate the "
        f"Model) writes `{METRICS_PATH.relative_to(REPO_ROOT)}`, this page "
        "will render AUROC/AUPRC curves, a confusion matrix, and a "
        "model-comparison table automatically."
    )
    st.stop()

with open(METRICS_PATH) as f:
    models = json.load(f)["models"]

st.subheader("Model comparison")
st.dataframe(
    [
        {
            "Model": name,
            "Accuracy": m["accuracy"],
            "Precision": m["precision"],
            "Recall": m["recall"],
            "F1": m["f1"],
            "AUROC": m["auroc"],
            "AUPRC": m["auprc"],
        }
        for name, m in models.items()
    ]
)

selected = st.selectbox("Model", list(models.keys()))
metrics = models[selected]

col1, col2 = st.columns(2)

with col1:
    st.subheader("ROC curve")
    roc_fig = go.Figure()
    roc_fig.add_trace(
        go.Scatter(x=metrics["fpr"], y=metrics["tpr"], mode="lines")
    )
    roc_fig.add_trace(
        go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(dash="dash"))
    )
    roc_fig.update_layout(
        xaxis_title="False Positive Rate", yaxis_title="True Positive Rate"
    )
    st.plotly_chart(roc_fig)

with col2:
    st.subheader("Precision-Recall curve")
    pr_fig = go.Figure()
    pr_fig.add_trace(
        go.Scatter(
            x=metrics["recall_curve"],
            y=metrics["precision_curve"],
            mode="lines",
        )
    )
    pr_fig.update_layout(xaxis_title="Recall", yaxis_title="Precision")
    st.plotly_chart(pr_fig)

st.subheader("Confusion matrix")
cm = metrics["confusion_matrix"]
cm_fig = go.Figure(
    go.Heatmap(
        z=cm,
        x=["Predicted negative", "Predicted positive"],
        y=["Actual negative", "Actual positive"],
        text=cm,
        texttemplate="%{text}",
    )
)
st.plotly_chart(cm_fig)
