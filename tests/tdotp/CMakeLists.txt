scene_add_test(TDoTP test_golden.py
  FILES tdotp_golden.cdf tdotp.dat
  ARGUMENTS -s ${CMAKE_BINARY_DIR}/bin/scene -i tdotp.dat -n tdotp.cdf -o tdotp_golden.cdf
  )
