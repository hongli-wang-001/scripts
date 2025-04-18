import argparse
import numpy as np
import utils
import model


def main():
    parser = argparse.ArgumentParser(description='Script to QC PurpleAir data, pair with AirNow data, and model.')
    parser.add_argument('--pa_dir', required=True, dest='pa_dir', type=str,
                        help='Directory for all the PurpleAir data. Should include sensor list.')
    parser.add_argument('--an_dir', required=True, dest='an_dir', type=str, help='Directory for all the AirNow data.')
    parser.add_argument('--output_dir', required=True, dest='output_dir', type=str,
                        help='Parent directory for all outputs.')
    parser.add_argument('--p', type=float, nargs='?', default=0.61,
                        help='Percentage difference threshold. Defaults to 0.61.')
    parser.add_argument('--a', type=float, nargs='?', default=5.,
                        help='Absolute difference threshold. Defaults to 5.')
    parser.add_argument('--m25', type=np.float32, nargs='?', default=3000,
                        help='Maximum threshold for PM2.5. Defaults to 3000.')
    parser.add_argument('--m10', type=np.float32, nargs='?', default=np.inf,
                        help='Maximum threshold for PM2.5. Defaults to infinity.')
    parser.add_argument('--r', type=float, nargs='?', default=500,
                        help='Collocation radius. Defaults to 500m.')
    parser.add_argument('--g', type=int, nargs='?', default=1,
                        help='0 if should skip generating hourly PurpleAir data.')
    args = parser.parse_args()
    if args.g:
        utils.generate_pa_hourly(args.pa_dir, args.output_dir, args.p, args.a, args.m25, args.m10)
    df = utils.pair_data(f"{args.output_dir}/hourly_qa_pm25", args.an_dir, args.r, "+proj=aea +lat_1=29.5 +lat_2=42.5")
    m = model.LinearRegression("PM25", ['pm25_pa', 'pm25_pa_2', 'humidity_a'])
    edges = np.concatenate(
        [[-1], np.arange(5, 101, 5), np.arange(110, 201, 10), np.arange(250, 401, 50), [df.PM25.max()]])
    edges = edges.reshape((-1, 1))
    m.fit(df, args.output_dir, edges)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
