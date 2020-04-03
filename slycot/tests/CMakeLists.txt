set(PYSOURCE

  __init__.py
  test.py
  test_ag08bd.py
  test_sb10jd.py
  test_sg02ad.py
  test_sg03ad.py
  test_tb05ad.py
  test_td04ad.py
  test_tg01ad.py
  test_tg01fd.py )

install(FILES ${PYSOURCE} DESTINATION slycot/tests)
