import logger

# rich_logger = logger.RichLogger(None)
# rich_logger.add_to_log("Hello,  World!", level = "info")
# rich_logger.add_to_log("Hello,  World!", level = "debug")
# rich_logger.add_to_log("Hello,  World!", level = "warning")
# rich_logger.add_to_log("Hello,  World!", level = "error")
# rich_logger.add_to_log("Hello,  World!", level = "critical")

simple_logger = logger.SimpleLogger("/Users/valehvpa/Desktop/test")
simple_logger.add_to_log("Hello, World!", level = "info")
simple_logger.add_to_log("Hello, World!", level = "debug")
simple_logger.add_to_log("Hello, World!", level = "warning")
simple_logger.add_to_log("Hello, World!", level = "error")
simple_logger.add_to_log("Hello, World!", level = "critical")