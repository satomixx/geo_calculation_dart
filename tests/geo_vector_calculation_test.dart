import 'package:test/test.dart';

import '../geo_vector_calculation.dart';

void main() {
  setUp(() {});

  group('vincentyDirect', () {
    test('valid', () {
      const lat = 24.288472; // 南鳥島の緯度
      const lon = 153.9707894; // 南鳥島の経度

      final result = vincentyDirect(
          latitude: lat,
          longitude: lon,
          azimuth: 276.8697566783211,
          distance: 3143772,
          ellipsoid: 1);
      if (result != null) {
        final _lat = result['latitude'];
        final _lon = result['longitude'];
        final _azimuth = result['azimuth'];
        print('latitude: $_lat');
        print('longitude: $_lon');
        print('azimuth: $_azimuth');

        expect(_lat, 24.455028332685178);
        expect(_lon, 122.90976683889555);
        expect(_azimuth, 83.78447885528264);
      }
    });
  });
}
