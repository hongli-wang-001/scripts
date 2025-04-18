from pathlib import Path
from typing import List

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn import linear_model


class LinearRegression:
    def __init__(self, target: str, covariates: List[str]):
        self.target = target
        self.covariates = covariates
        self.model = None

    def fit(self, df: pd.DataFrame, save_dir: str, edges: np.ndarray):
        print("Fitting linear regression.")
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        d = df.groupby(['AQSID', 'datetime_utc'])[[self.target] + self.covariates].agg(np.nanmean)
        d = d.dropna()
        x = d[self.covariates].to_numpy()
        y = d[[self.target]].to_numpy()

        self.model = linear_model.LinearRegression()
        self.model.fit(x, y)
        qa_pred = self.model.predict(x)
        for var, coef in zip(self.covariates, self.model.coef_[0]):
            print(f"{var}: {coef}")
        print(f"intercept: {self.model.intercept_}")
        print(f"R2: {self.model.score(x, y)}")
        print(f"RMSE: {np.mean((y - qa_pred) ** 2) ** 0.5}")

        print(f"Saving outputs to {save_dir}.")

        df['pred'] = self.model.predict(df[self.covariates].to_numpy())
        df['residuals'] = df.pred - df.PM25
        an = df.PM25.to_numpy()
        df['bin'] = np.argmax((an > edges[:-1]) & (an <= edges[1:]), axis=0)
        res_stats = df.groupby('bin').agg(res_mean=("residuals", "mean"),
                                          res_std=("residuals", "std"))
        edges[0, 0] = 0
        with open(save_dir / Path("residuals.txt"), 'w') as f:
            f.write("low, high, mean, std\n")
            f.writelines([f"{low}, {high}, {mean}, {std}\n"
                          for low, high, mean, std
                          in zip(edges[:-1, 0], edges[1:, 0], res_stats.res_mean, res_stats.res_std)])

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        fig.suptitle(f"PM2.5 Correlation (PurpleAir/Modeled vs. AirNow)")
        ax1.set_title("PurpleAir PM2.5 vs. AirNow PM2.5")
        sns.kdeplot(x=d.pm25_pa, y=y.flatten(), cmap="terrain", fill=True, bw_adjust=.5, thresh=0.01, clip=(0, 400),
                    ax=ax1)
        ax1.set_xlabel("PurpleAir PM2.5")
        ax1.set_ylabel("AirNow PM2.5")
        ax1.plot([0, 400], [0, 400], color='r')
        ax2.set_title("Predicted PM2.5 vs AirNow PM2.5")
        sns.kdeplot(x=qa_pred.flatten(), y=y.flatten(), cmap="terrain", fill=True, bw_adjust=.5, thresh=0.01,
                    clip=(0, 400),
                    ax=ax2)
        ax2.set_xlabel("Predicted PM2.5")
        ax2.set_ylabel("AirNow PM2.5")
        ax2.plot([0, 400], [0, 400], color='r')
        plt.savefig(save_dir / Path("fit.png"))
        plt.close()

        with open(save_dir / Path("lin_reg.txt"), 'w') as f:
            f.writelines([f'{var}, {coef}\n' for var, coef in zip(self.covariates, self.model.coef_[0])])

        station_dir = save_dir / Path("images/") / Path("station_plots")
        station_dir.mkdir(parents=True, exist_ok=True)
        index = pd.date_range(df.datetime_utc.min(), df.datetime_utc.max(), freq="h", name="datetime_utc")
        for sid, data in df.groupby('AQSID'):
            graph_params = data.groupby("datetime_utc").agg(an_pm=("PM25", "first"),
                                                            pred_mean=("pred", "mean"),
                                                            pred_low=("pred", "min"),
                                                            pred_high=("pred", "max")).reindex(index)
            plt.figure(figsize=(12, 8))
            plt.xticks(rotation=45, ha='right')
            plt.title(f"({sid}) {data.SiteName.iloc[0]} {data.StateName.iloc[0]}, {data.CountryCode.iloc[0]}")
            plt.ylabel("PM2.5 concentration (µg/m3)")
            plt.ylabel("time")
            plt.plot(graph_params.index, graph_params.pred_mean, lw=0.75)
            plt.fill_between(x=graph_params.index, y1=graph_params.pred_low, y2=graph_params.pred_high)
            plt.semilogy(graph_params.index, graph_params.an_pm, color="r", lw=0.75)
            plt.savefig(station_dir / Path(f"{sid}.png"))
            plt.close()
