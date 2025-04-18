import glob
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd
import pytz
import timezonefinder
from matplotlib import pyplot as plt
from pandas.errors import ParserError

tf = timezonefinder.TimezoneFinder()


def get_utc(x):
    tz = pytz.timezone(tf.timezone_at(lat=x.latitude, lng=x.longitude))
    return x.datetime_local - tz.localize(x.datetime_local).utcoffset()


def generate_pa_hourly(data_dir: str,
                       save_dir: str,
                       percentage_diff_threshold: float,
                       abs_diff_threshold: float,
                       pm25_max: float,
                       pm10_max: float) -> None:
    print(f"Generating PurpleAir hourly data at {save_dir}.")
    sensor_list = pd.read_csv(f"{data_dir}/sensors_list.csv", quotechar='"').drop("name", axis=1)
    sensor_list = sensor_list.set_axis(['sensor_index', 'location_type', 'latitude', 'longitude'], axis=1)
    pa_all = None
    for fn in glob.glob(f"{data_dir}/sensorsID*.csv"):
        curr = pd.read_csv(fn)
        if pa_all is None:
            pa_all = curr
        else:
            pa_all = pd.concat([pa_all, curr])
    pa_all = pd.merge(pa_all, sensor_list, on='sensor_index', how='inner')
    pa_all['datetime_local'] = pd.to_datetime(pa_all['date_time_utc'], format="%Y-%m-%d %H:%M:%S")
    pa_all = pa_all.drop(["date_time_utc"], axis=1)
    pa_all['datetime_utc'] = pa_all.apply(get_utc, axis=1)

    pa_all['pm25_pa'] = (pa_all['pm2.5_cf_1_a'] + pa_all['pm2.5_cf_1_b']) / 2
    pa_all['pm10_pa'] = (pa_all['pm10.0_cf_1_a'] + pa_all['pm10.0_cf_1_b']) / 2

    pa_all_qa_pm25 = pa_all.copy()
    pa_all_qa_pm25 = pa_all_qa_pm25[
        (pa_all_qa_pm25['pm2.5_cf_1_a'] <= pm25_max) & (pa_all_qa_pm25['pm2.5_cf_1_b'] <= pm25_max)]
    pa_all_qa_pm25 = pa_all_qa_pm25[(pa_all_qa_pm25['pm2.5_cf_1_a'] >= 0) & (pa_all_qa_pm25['pm2.5_cf_1_b'] >= 0)]
    pa_all_qa_pm25 = pa_all_qa_pm25[pa_all_qa_pm25["pm2.5_cf_1_a"] != pa_all_qa_pm25["pm2.5_cf_1_b"]]
    pa_all_qa_pm25 = pa_all_qa_pm25[pa_all_qa_pm25["pm2.5_cf_1_a"] != pa_all_qa_pm25["pm10.0_cf_1_a"]]
    pa_all_qa_pm25 = pa_all_qa_pm25[pa_all_qa_pm25["pm2.5_cf_1_b"] != pa_all_qa_pm25["pm10.0_cf_1_b"]]
    pa_all_qa_pm25 = pa_all_qa_pm25[pa_all_qa_pm25["temperature_a"] < 1000]
    pa_all_qa_pm25 = pa_all_qa_pm25[(pa_all_qa_pm25["humidity_a"] >= 0) & (pa_all_qa_pm25["humidity_a"] <= 100)]
    abs_diff_25 = (pa_all_qa_pm25["pm2.5_cf_1_a"] - pa_all_qa_pm25["pm2.5_cf_1_b"]).abs()
    per_diff_25 = abs_diff_25 * 2 / (pa_all_qa_pm25["pm2.5_cf_1_a"] + pa_all_qa_pm25["pm2.5_cf_1_b"])
    pa_all_qa_pm25 = pa_all_qa_pm25[(abs_diff_25 < abs_diff_threshold) | (per_diff_25 < percentage_diff_threshold)]

    pa_all_qa_pm10 = pa_all.copy()
    pa_all_qa_pm10 = pa_all_qa_pm10[
        (pa_all_qa_pm10['pm10.0_cf_1_a'] <= pm10_max) & (pa_all_qa_pm10['pm10.0_cf_1_b'] <= pm10_max)]
    pa_all_qa_pm10 = pa_all_qa_pm10[(pa_all_qa_pm10['pm10.0_cf_1_a'] >= 0) & (pa_all_qa_pm10['pm10.0_cf_1_b'] >= 0)]
    pa_all_qa_pm10 = pa_all_qa_pm10[pa_all_qa_pm10["pm10.0_cf_1_a"] != pa_all_qa_pm10["pm10.0_cf_1_b"]]
    pa_all_qa_pm10 = pa_all_qa_pm10[pa_all_qa_pm10["pm2.5_cf_1_a"] != pa_all_qa_pm10["pm10.0_cf_1_a"]]
    pa_all_qa_pm10 = pa_all_qa_pm10[pa_all_qa_pm10["pm2.5_cf_1_b"] != pa_all_qa_pm10["pm10.0_cf_1_b"]]
    pa_all_qa_pm10 = pa_all_qa_pm10[pa_all_qa_pm10["temperature_a"] < 1000]
    pa_all_qa_pm10 = pa_all_qa_pm10[(pa_all_qa_pm10["humidity_a"] >= 0) & (pa_all_qa_pm10["humidity_a"] <= 100)]
    abs_diff_10 = (pa_all_qa_pm10["pm10.0_cf_1_a"] - pa_all_qa_pm10["pm10.0_cf_1_b"]).abs()
    per_diff_10 = abs_diff_10 * 2 / (pa_all_qa_pm10["pm10.0_cf_1_a"] + pa_all_qa_pm10["pm10.0_cf_1_b"])
    pa_all_qa_pm10 = pa_all_qa_pm10[(abs_diff_10 < abs_diff_threshold) | (per_diff_10 < percentage_diff_threshold)]

    # pa_all.to_csv(f"{save_dir}/all.csv", index=False, date_format="%Y%m%d %H:%M:%S")
    # pa_all_qa_pm25.to_csv(f"{save_dir}/all_qa_pm25.csv", index=False, date_format="%Y%m%d %H:%M:%S")
    # pa_all_qa_pm10.to_csv(f"{save_dir}/all_qa_pm10.csv", index=False, date_format="%Y%m%d %H:%M:%S")

    (Path(save_dir) / Path("hourly")).mkdir(parents=True, exist_ok=True)
    (Path(save_dir) / Path("hourly_qa_pm25")).mkdir(parents=True, exist_ok=True)
    (Path(save_dir) / Path("hourly_qa_pm10")).mkdir(parents=True, exist_ok=True)
    for time, group in pa_all.groupby("datetime_utc"):
        hourly = group.sort_values("sensor_index").reset_index(drop=True)
        hourly.to_csv(f"{save_dir}/hourly/{time.strftime('%Y%m%d_%H')}.csv", index=False,
                      date_format="%Y%m%d %H:%M:%S")

    for time, group in pa_all_qa_pm25.groupby("datetime_utc"):
        hourly = group.sort_values("sensor_index").reset_index(drop=True)
        hourly.to_csv(f"{save_dir}/hourly_qa_pm25/{time.strftime('%Y%m%d_%H')}.csv", index=False,
                      date_format="%Y%m%d %H:%M:%S")

    for time, group in pa_all_qa_pm10.groupby("datetime_utc"):
        hourly = group.sort_values("sensor_index").reset_index(drop=True)
        hourly.to_csv(f"{save_dir}/hourly_qa_pm10/{time.strftime('%Y%m%d_%H')}.csv", index=False,
                      date_format="%Y%m%d %H:%M:%S")

    fig = plt.figure(constrained_layout=True)
    fig.suptitle("PurpleAir PM histograms")
    subfigs = fig.subfigures(nrows=2, ncols=1)
    subfig_titles = ["PM2.5", "PM10"]
    subplot_titles = ["All", "QC"]
    subfig_data = [[pa_all.pm25_pa, pa_all_qa_pm25.pm25_pa], [pa_all.pm10_pa, pa_all_qa_pm10.pm10_pa]]
    for row, subfig in enumerate(subfigs):
        subfig.suptitle(subfig_titles[row])
        axs = subfig.subplots(nrows=1, ncols=2)
        for col, ax in enumerate(axs):
            ax.hist(subfig_data[row][col], bins=20, log=True)
            ax.set_title(subplot_titles[col])
    image_path = Path(save_dir) / Path("images")
    image_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(image_path / Path("pa_hist.png"))
    plt.close()

    fns = sorted(glob.glob(f"{save_dir}/hourly_qa_pm25/*.csv"))
    start_date = pd.to_datetime(fns[0][-15:-4], format="%Y%m%d_%H")
    end_date = pd.to_datetime(fns[-1][-15:-4], format="%Y%m%d_%H")

    fig, axs = plt.subplots(2, 2, figsize=(16, 16))
    dates_pm25 = pd.date_range(start_date, end_date, freq="H")
    all_count_pm25 = []
    qa_count_pm25 = []
    for date in dates_pm25:
        try:
            all = pd.read_csv(f"{save_dir}/hourly/{date.strftime('%Y%m%d_%H')}.csv")
            all_count_pm25.append(len(all))
        except FileNotFoundError:
            all_count_pm25.append(0)
        try:
            qa = pd.read_csv(f"{save_dir}/hourly_qa_pm25/{date.strftime('%Y%m%d_%H')}.csv")
            qa_count_pm25.append(len(qa))
        except FileNotFoundError:
            qa_count_pm25.append(0)
    axs[0, 0].plot(dates_pm25, all_count_pm25, label="all")
    axs[0, 0].plot(dates_pm25, qa_count_pm25, label="qa")
    axs[0, 0].set_title("Number of Data Records (PM2.5)")
    axs[0, 0].legend()
    axs[0, 0].set_xlabel("date")
    axs[0, 0].set_xticks(axs[0, 0].get_xticks(), axs[0, 0].get_xticklabels(), rotation=45, ha='right')
    axs[0, 0].set_ylabel("count")

    axs[0, 1].plot(dates_pm25, [q / a for q, a in zip(qa_count_pm25, all_count_pm25)])
    axs[0, 1].set_title("Ratio of QA Record Count to Total Record Count (PM2.5)")
    axs[0, 1].set_xlabel("date")
    axs[0, 1].set_ylabel("ratio")
    axs[0, 1].set_xticks(axs[0, 1].get_xticks(), axs[0, 1].get_xticklabels(), rotation=45, ha='right')

    dates_pm10 = pd.date_range(start_date, end_date, freq="H")
    all_count_pm10 = []
    qa_count_pm10 = []
    for date in dates_pm10:
        try:
            all = pd.read_csv(f"{save_dir}/hourly/{date.strftime('%Y%m%d_%H')}.csv")
            all_count_pm10.append(len(all))
        except FileNotFoundError:
            all_count_pm10.append(0)
        try:
            qa = pd.read_csv(f"{save_dir}/hourly_qa_pm10/{date.strftime('%Y%m%d_%H')}.csv")
            qa_count_pm10.append(len(qa))
        except FileNotFoundError:
            qa_count_pm10.append(0)
    axs[1, 0].plot(dates_pm10, all_count_pm10, label="all")
    axs[1, 0].plot(dates_pm10, qa_count_pm10, label="qa")
    axs[1, 0].set_title("Number of Data Records (PM10)")
    axs[1, 0].legend()
    axs[1, 0].set_xlabel("date")
    axs[1, 0].set_xticks(axs[1, 0].get_xticks(), axs[1, 0].get_xticklabels(), rotation=45, ha='right')
    axs[1, 0].set_ylabel("count")

    axs[1, 1].plot(dates_pm10, [q / a for q, a in zip(qa_count_pm10, all_count_pm10)])
    axs[1, 1].set_title("Ratio of QA Record Count to Total Record Count (PM10)")
    axs[1, 1].set_xlabel("date")
    axs[1, 1].set_ylabel("ratio")
    axs[1, 1].set_xticks(axs[1, 1].get_xticks(), axs[1, 1].get_xticklabels(), rotation=45, ha='right')

    plt.savefig(image_path / Path("record_count.png"))
    plt.close()

    (image_path / Path("qa_maps")).mkdir(parents=True, exist_ok=True)
    for date in dates_pm25:
        plt.figure(figsize=(12, 10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.STATES, linewidth=0.3)
        d = pa_all[pa_all.datetime_utc == date]
        dqa = pa_all_qa_pm25[pa_all_qa_pm25.datetime_utc == date]
        plt.scatter(d.longitude, d.latitude, transform=ccrs.PlateCarree(), color="r", s=0.5)
        plt.scatter(dqa.longitude, dqa.latitude, transform=ccrs.PlateCarree(), color="b", s=0.5)
        plt.title(date.strftime('%Y-%m-%d- %H:%M:%S'))
        plt.savefig(image_path / Path("qa_maps") / Path(date.strftime('%Y%m%d%H.png')))
        plt.close()


def pair_data(pa_hourly_dir: str,
              an_hourly_dir: str,
              collocation_radius: float,
              crs: str) -> pd.DataFrame:
    print(f"Collocating data ({collocation_radius}m).")
    start_date = min(pd.to_datetime(sorted(os.listdir(pa_hourly_dir))[0], format="%Y%m%d_%H.csv"),
                     pd.to_datetime(sorted(os.listdir(an_hourly_dir))[0], format="HourlyAQObs_%Y%m%d%H.dat"))
    end_date = max(pd.to_datetime(sorted(os.listdir(pa_hourly_dir))[-1], format="%Y%m%d_%H.csv"),
                     pd.to_datetime(sorted(os.listdir(an_hourly_dir))[-1], format="HourlyAQObs_%Y%m%d%H.dat"))
    date_range = pd.date_range(start=start_date, end=end_date, freq="H")

    # all = [None] * len(collocation_radius)
    all_qa = None
    for date in date_range:
        try:
            # pa = pd.read_csv(f"{pa_hourly_dir}hourly/{date.strftime('%Y%m%d_%H')}.csv")
            pa_qa = pd.read_csv(f"{pa_hourly_dir}/{date.strftime('%Y%m%d_%H')}.csv")
            an = pd.read_csv(f"{an_hourly_dir}/HourlyAQObs_{date.strftime('%Y%m%d%H')}.dat")
            an = an.loc[(an.Latitude < 52) & (an.Latitude > 24) & (an.Longitude < -64) & (an.Longitude > -137) & (
                ~np.isnan(an.PM25)) & (an.PM25 >= 0)]
            # pa = gpd.GeoDataFrame(data=pa, geometry=gpd.points_from_xy(pa.longitude, pa.latitude),
            #                       crs="EPSG:4326").to_crs(crs)
            pa_qa = gpd.GeoDataFrame(data=pa_qa, geometry=gpd.points_from_xy(pa_qa.longitude, pa_qa.latitude),
                                     crs="EPSG:4326").to_crs(crs)
            an = gpd.GeoDataFrame(data=an, geometry=gpd.points_from_xy(an.Longitude, an.Latitude),
                                  crs="EPSG:4326").to_crs(crs)
            an.geometry = an.geometry.buffer(collocation_radius)
            # joined = gpd.sjoin_nearest(pa, an, how="inner", max_distance=0.001)
            # joined = joined.drop("index_right", axis=1)
            joined_qa = gpd.sjoin_nearest(pa_qa, an, how="inner", max_distance=0.001)
            joined_qa = joined_qa.drop("index_right", axis=1)
            # if all[i] is None:
            #     all[i] = joined
            # else:
            #     all[i] = pd.concat([all[i], joined])
            if all_qa is None:
                all_qa = joined_qa
            else:
                all_qa = pd.concat([all_qa, joined_qa])
        except (FileNotFoundError, ParserError) as e:
            print(f"failed to process {date} {e}")
    return all_qa
