set(PYSOURCE

  __init__.py test.py test_sg02ad.py test_sg03ad.py test_tb05ad.py
  test_td04ad.py)

install(FILES ${PYSOURCE} DESTINATION slycot/tests)
