set(PYSOURCE

  __init__.py
  test_ab08n.py
  test_ag08bd.py
  test_examples.py
  test_mb.py
  test_mc.py
  test_sb.py
  test_sg02ad.py
  test_sg03ad.py
  test_tb05ad.py
  test_td04ad.py
  test_tg01ad.py
  test_tg01fd.py )

install(FILES ${PYSOURCE}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE
                    WORLD_READ WORLD_EXECUTE
        DESTINATION slycot/tests)
