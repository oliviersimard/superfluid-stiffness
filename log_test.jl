using Logging

#default:
#Logging.configure(level=WARNING)

Logging.configure(level=DEBUG)

function log_test()
    debug("debug message")
    info("info message")
    warn("warning message")
    err("error message")
    critical("critical message")
end


