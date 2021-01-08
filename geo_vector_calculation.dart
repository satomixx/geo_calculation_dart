import 'dart:math';
import 'package:vector_math/vector_math.dart';

// ref: http://hhsprings.pinoko.jp/personal_works/site-hhs_downloads/latitude_ja.pdf

/// 楕円体
const ellispoidGrs80 = 1; // GRS80
const ellispoidWgs84 = 2; // WGS84

/// 楕円体別の長軸半径と扁平率
const geodeticDataEllispoidGrs80 = [
  6378137.0, // [GRS80]長軸半径
  1 / 298.257222101, // [GRS80]扁平率
];
const geodeticDataEllispoidWgs84 = [
  6378137.0, // [WGS84]長軸半径
  1 / 298.257223563, // [WGS84]扁平率
];

// 反復計算の上限回数
const iterationLimit = 1000;

/// Vincenty法(順解法)
/// 始点の座標(緯度経度)と方位角と距離から、終点の座標と方位角を求める
/// latitude 緯度
/// longitude 経度
/// azimuth 方位角
/// distance 距離
/// ellipsoid 楕円体
/// return 終点の座標、方位角
Map<String, double> vincentyDirect(
    {double latitude,
    double longitude,
    double azimuth,
    double distance,
    int ellipsoid}) {
  // 計算時に必要な長軸半径(majorAxisRadius)と扁平率(flatness)を定数から取得し、
  // 短軸半径(minorAxisRadius)を算出する
  // 楕円体が未指定の場合はGRS80の値を用いる
  var majorAxisRadius = geodeticDataEllispoidGrs80[0];
  var flatness = geodeticDataEllispoidGrs80[1];
  if (ellipsoid != null && ellipsoid == 2) {
    majorAxisRadius = geodeticDataEllispoidWgs84[0];
    flatness = geodeticDataEllispoidWgs84[0];
  }

  final minorAxisRadius = (1 - flatness) * majorAxisRadius;

  // ラジアンに変換する(距離以外)
  final phi1 = radians(latitude); // φ, phi
  final lambda1 = radians(longitude); // λ, lambda
  final alpha1 = radians(azimuth); // α, alpha
  // final s = distance;

  final sinAlpha1 = sin(alpha1);
  final cosAlpha1 = cos(alpha1);

  // 更成緯度(補助球上の緯度)
  final reducedLatitude = atan((1 - flatness) * tan(phi1));

  final sinRl = sin(reducedLatitude);
  final cosRl = cos(reducedLatitude);
  final tanRl = tan(reducedLatitude);

  final sigma1 = atan2(tanRl, cosAlpha1); // σ. sigma
  final sinAlpha = cosRl * sinAlpha1;
  final cos2alpha = 1 - pow(sinAlpha, 2);

  final _p = pow(majorAxisRadius, 2) - pow(minorAxisRadius, 2);
  final u2 = cos2alpha * _p / pow(minorAxisRadius, 2);
  final A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
  final B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));

  // σをdistance/(minorAxisRadius*A)で初期化
  var sigma = distance / (minorAxisRadius * A); // σ, sigma

  // 以下の計算をσが収束するまで反復する
  // 地点によっては収束しないことがあり得るため、反復回数に上限を設ける
  final cos2SigmaM = cos(2 * sigma1 + sigma);
  final sinSigma = sin(sigma);
  final cosSigma = cos(sigma);
  for (var i = 0; i < iterationLimit; i++) {
    final lambdaSigma = B *
        sinSigma *
        (cos2SigmaM +
            B /
                4 *
                (cosSigma * (-1 + 2 * pow(cos2SigmaM, 2)) -
                    B /
                        6 *
                        cos2SigmaM *
                        (-3 + 4 * pow(sinSigma, 2)) *
                        (-3 + 4 * pow(cos2SigmaM, 2))));
    final sigmaDash = sigma;
    sigma = distance / (minorAxisRadius * A) + lambdaSigma;

    // 偏差が.000000000001以下ならbreak
    if ((sigma - sigmaDash).abs() <= 1e-12) {
      break;
    }
    // else:
    //     // 計算が収束しなかった場合はnullを返す
    //     return null;
  }

  // σが所望の精度まで収束したら以下の計算を行う
  final x = sinRl * sinSigma - cosRl * cosSigma * cosAlpha1;
  final phi2 = atan2(sinRl * cosSigma + cosRl * sinSigma * cosAlpha1,
      (1 - flatness) * sqrt(pow(sinAlpha, 2) + pow(x, 2)));
  final lambda = atan2(
      sinSigma * sinAlpha1, cosRl * cosSigma - sinRl * sinSigma * cosAlpha1);
  final C = flatness / 16 * cos2alpha * (4 + flatness * (4 - 3 * cos2alpha));
  final L = lambda -
      (1 - C) *
          flatness *
          sinAlpha *
          (sigma +
              C *
                  sinSigma *
                  (cos2SigmaM + C * cosSigma * (-1 + 2 * pow(cos2SigmaM, 2))));
  final lambda2 = L + lambda1;

  final alpha2 = atan2(sinAlpha, -x) + pi;

  return {
    'latitude': degrees(phi2),
    'longitude': degrees(lambda2),
    'azimuth': degrees(alpha2),
  };
}

void main(List<String> args) {
  args.forEach(print);
  final _latitude = double.parse(args[1]);
  final _longitude = double.parse(args[2]);
  final _azimuth = double.parse(args[3]);
  final _distance = double.parse(args[4]);

  if (args[5] == null) {
    vincentyDirect(
      latitude: _latitude,
      longitude: _longitude,
      azimuth: _azimuth,
      distance: _distance,
    );
  } else {
    vincentyDirect(
        latitude: _latitude,
        longitude: _longitude,
        azimuth: _azimuth,
        distance: _distance,
        ellipsoid: int.parse(args[5]));
  }
}
