#include <iostream>
#include <thread>
#include <vector>
#include <string>
#include <cstring>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <fcntl.h>
#include <errno.h>

class TunnelServer {
private:
    int listen_socket;
    std::string target_host;
    int target_port;
    int listen_port;
    bool is_ipv6_listen;
    bool is_ipv6_target;

public:
    TunnelServer(const std::string& listen_addr, int listen_port, 
                 const std::string& target_host, int target_port)
        : target_host(target_host), target_port(target_port), 
          listen_port(listen_port), listen_socket(-1) {
        
        // Determine if listen address is IPv6
        is_ipv6_listen = (listen_addr.find(':') != std::string::npos);
        
        // Determine if target is IPv6
        is_ipv6_target = (target_host.find(':') != std::string::npos);
        
        createListenSocket(listen_addr);
    }
    
    ~TunnelServer() {
        if (listen_socket >= 0) {
            close(listen_socket);
        }
    }

private:
    void createListenSocket(const std::string& listen_addr) {
        if (is_ipv6_listen) {
            listen_socket = socket(AF_INET6, SOCK_STREAM, 0);
            if (listen_socket < 0) {
                throw std::runtime_error("Failed to create IPv6 socket");
            }
            
            // Allow IPv4 connections on IPv6 socket
            int no = 0;
            setsockopt(listen_socket, IPPROTO_IPV6, IPV6_V6ONLY, &no, sizeof(no));
            
            struct sockaddr_in6 addr;
            memset(&addr, 0, sizeof(addr));
            addr.sin6_family = AF_INET6;
            addr.sin6_port = htons(listen_port);
            
            if (inet_pton(AF_INET6, listen_addr.c_str(), &addr.sin6_addr) <= 0) {
                throw std::runtime_error("Invalid IPv6 address");
            }
            
            int reuse = 1;
            setsockopt(listen_socket, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse));
            
            if (bind(listen_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
                throw std::runtime_error("Failed to bind IPv6 socket");
            }
        } else {
            listen_socket = socket(AF_INET, SOCK_STREAM, 0);
            if (listen_socket < 0) {
                throw std::runtime_error("Failed to create IPv4 socket");
            }
            
            struct sockaddr_in addr;
            memset(&addr, 0, sizeof(addr));
            addr.sin_family = AF_INET;
            addr.sin_port = htons(listen_port);
            
            if (inet_pton(AF_INET, listen_addr.c_str(), &addr.sin_addr) <= 0) {
                throw std::runtime_error("Invalid IPv4 address");
            }
            
            int reuse = 1;
            setsockopt(listen_socket, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse));
            
            if (bind(listen_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
                throw std::runtime_error("Failed to bind IPv4 socket");
            }
        }
        
        if (listen(listen_socket, 10) < 0) {
            throw std::runtime_error("Failed to listen on socket");
        }
    }
    
    int connectToTarget() {
        int target_socket;
        
        if (is_ipv6_target) {
            target_socket = socket(AF_INET6, SOCK_STREAM, 0);
            if (target_socket < 0) {
                return -1;
            }
            
            struct sockaddr_in6 addr;
            memset(&addr, 0, sizeof(addr));
            addr.sin6_family = AF_INET6;
            addr.sin6_port = htons(target_port);
            
            if (inet_pton(AF_INET6, target_host.c_str(), &addr.sin6_addr) <= 0) {
                close(target_socket);
                return -1;
            }
            
            if (connect(target_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
                close(target_socket);
                return -1;
            }
        } else {
            target_socket = socket(AF_INET, SOCK_STREAM, 0);
            if (target_socket < 0) {
                return -1;
            }
            
            struct sockaddr_in addr;
            memset(&addr, 0, sizeof(addr));
            addr.sin_family = AF_INET;
            addr.sin_port = htons(target_port);
            
            if (inet_pton(AF_INET, target_host.c_str(), &addr.sin_addr) <= 0) {
                close(target_socket);
                return -1;
            }
            
            if (connect(target_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
                close(target_socket);
                return -1;
            }
        }
        
        return target_socket;
    }
    
    void handleConnection(int client_socket) {
        int target_socket = connectToTarget();
        if (target_socket < 0) {
            std::cerr << "Failed to connect to target" << std::endl;
            close(client_socket);
            return;
        }
        
        std::cout << "Tunnel established" << std::endl;
        
        // Create threads for bidirectional data transfer
        std::thread t1([this, client_socket, target_socket]() {
            this->transferData(client_socket, target_socket, "Client->Target");
        });
        
        std::thread t2([this, client_socket, target_socket]() {
            this->transferData(target_socket, client_socket, "Target->Client");
        });
        
        t1.join();
        t2.join();
        
        close(client_socket);
        close(target_socket);
        std::cout << "Tunnel closed" << std::endl;
    }
    
    void transferData(int from_socket, int to_socket, const std::string& direction) {
        char buffer[4096];
        ssize_t bytes_read, bytes_written;
        
        while (true) {
            bytes_read = recv(from_socket, buffer, sizeof(buffer), 0);
            if (bytes_read <= 0) {
                if (bytes_read < 0) {
                    std::cerr << "Error reading from " << direction << ": " 
                              << strerror(errno) << std::endl;
                }
                break;
            }
            
            bytes_written = 0;
            while (bytes_written < bytes_read) {
                ssize_t result = send(to_socket, buffer + bytes_written, 
                                    bytes_read - bytes_written, MSG_NOSIGNAL);
                if (result <= 0) {
                    if (result < 0) {
                        std::cerr << "Error writing to " << direction << ": " 
                                  << strerror(errno) << std::endl;
                    }
                    return;
                }
                bytes_written += result;
            }
        }
    }

public:
    void start() {
        std::cout << "Tunnel server started on port " << listen_port << std::endl;
        std::cout << "Tunneling to " << target_host << ":" << target_port << std::endl;
        
        while (true) {
            struct sockaddr_storage client_addr;
            socklen_t client_len = sizeof(client_addr);
            
            int client_socket = accept(listen_socket, (struct sockaddr*)&client_addr, &client_len);
            if (client_socket < 0) {
                std::cerr << "Accept failed: " << strerror(errno) << std::endl;
                continue;
            }
            
            // Handle each connection in a separate thread
            std::thread([this, client_socket]() {
                this->handleConnection(client_socket);
            }).detach();
        }
    }
};

class TunnelClient {
private:
    std::string local_host;
    int local_port;
    std::string tunnel_host;
    int tunnel_port;
    bool is_ipv6_local;
    bool is_ipv6_tunnel;

public:
    TunnelClient(const std::string& local_host, int local_port,
                 const std::string& tunnel_host, int tunnel_port)
        : local_host(local_host), local_port(local_port),
          tunnel_host(tunnel_host), tunnel_port(tunnel_port) {
        
        is_ipv6_local = (local_host.find(':') != std::string::npos);
        is_ipv6_tunnel = (tunnel_host.find(':') != std::string::npos);
    }

private:
    int connectToTunnel() {
        int tunnel_socket;
        
        if (is_ipv6_tunnel) {
            tunnel_socket = socket(AF_INET6, SOCK_STREAM, 0);
            if (tunnel_socket < 0) {
                return -1;
            }
            
            struct sockaddr_in6 addr;
            memset(&addr, 0, sizeof(addr));
            addr.sin6_family = AF_INET6;
            addr.sin6_port = htons(tunnel_port);
            
            if (inet_pton(AF_INET6, tunnel_host.c_str(), &addr.sin6_addr) <= 0) {
                close(tunnel_socket);
                return -1;
            }
            
            if (connect(tunnel_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
                close(tunnel_socket);
                return -1;
            }
        } else {
            tunnel_socket = socket(AF_INET, SOCK_STREAM, 0);
            if (tunnel_socket < 0) {
                return -1;
            }
            
            struct sockaddr_in addr;
            memset(&addr, 0, sizeof(addr));
            addr.sin_family = AF_INET;
            addr.sin_port = htons(tunnel_port);
            
            if (inet_pton(AF_INET, tunnel_host.c_str(), &addr.sin_addr) <= 0) {
                close(tunnel_socket);
                return -1;
            }
            
            if (connect(tunnel_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
                close(tunnel_socket);
                return -1;
            }
        }
        
        return tunnel_socket;
    }
    
    int bindLocalSocket() {
        int local_socket;
        
        if (is_ipv6_local) {
            local_socket = socket(AF_INET6, SOCK_STREAM, 0);
            if (local_socket < 0) {
                return -1;
            }
            
            // Allow IPv4 connections on IPv6 socket
            int no = 0;
            setsockopt(local_socket, IPPROTO_IPV6, IPV6_V6ONLY, &no, sizeof(no));
            
            struct sockaddr_in6 addr;
            memset(&addr, 0, sizeof(addr));
            addr.sin6_family = AF_INET6;
            addr.sin6_port = htons(local_port);
            
            if (inet_pton(AF_INET6, local_host.c_str(), &addr.sin6_addr) <= 0) {
                close(local_socket);
                return -1;
            }
            
            int reuse = 1;
            setsockopt(local_socket, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse));
            
            if (bind(local_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
                close(local_socket);
                return -1;
            }
        } else {
            local_socket = socket(AF_INET, SOCK_STREAM, 0);
            if (local_socket < 0) {
                return -1;
            }
            
            struct sockaddr_in addr;
            memset(&addr, 0, sizeof(addr));
            addr.sin_family = AF_INET;
            addr.sin_port = htons(local_port);
            
            if (inet_pton(AF_INET, local_host.c_str(), &addr.sin_addr) <= 0) {
                close(local_socket);
                return -1;
            }
            
            int reuse = 1;
            setsockopt(local_socket, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse));
            
            if (bind(local_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
                close(local_socket);
                return -1;
            }
        }
        
        if (listen(local_socket, 10) < 0) {
            close(local_socket);
            return -1;
        }
        
        return local_socket;
    }
    
    void handleConnection(int local_client_socket) {
        int tunnel_socket = connectToTunnel();
        if (tunnel_socket < 0) {
            std::cerr << "Failed to connect to tunnel server" << std::endl;
            close(local_client_socket);
            return;
        }
        
        std::cout << "Tunnel connection established" << std::endl;
        
        // Create threads for bidirectional data transfer
        std::thread t1([this, local_client_socket, tunnel_socket]() {
            this->transferData(local_client_socket, tunnel_socket, "Local->Tunnel");
        });
        
        std::thread t2([this, local_client_socket, tunnel_socket]() {
            this->transferData(tunnel_socket, local_client_socket, "Tunnel->Local");
        });
        
        t1.join();
        t2.join();
        
        close(local_client_socket);
        close(tunnel_socket);
        std::cout << "Tunnel connection closed" << std::endl;
    }
    
    void transferData(int from_socket, int to_socket, const std::string& direction) {
        char buffer[4096];
        ssize_t bytes_read, bytes_written;
        
        while (true) {
            bytes_read = recv(from_socket, buffer, sizeof(buffer), 0);
            if (bytes_read <= 0) {
                if (bytes_read < 0) {
                    std::cerr << "Error reading from " << direction << ": " 
                              << strerror(errno) << std::endl;
                }
                break;
            }
            
            bytes_written = 0;
            while (bytes_written < bytes_read) {
                ssize_t result = send(to_socket, buffer + bytes_written, 
                                    bytes_read - bytes_written, MSG_NOSIGNAL);
                if (result <= 0) {
                    if (result < 0) {
                        std::cerr << "Error writing to " << direction << ": " 
                                  << strerror(errno) << std::endl;
                    }
                    return;
                }
                bytes_written += result;
            }
        }
    }

public:
    void start() {
        int local_socket = bindLocalSocket();
        if (local_socket < 0) {
            throw std::runtime_error("Failed to create local listening socket");
        }
        
        std::cout << "Tunnel client started on " << local_host << ":" << local_port << std::endl;
        std::cout << "Connecting through tunnel at " << tunnel_host << ":" << tunnel_port << std::endl;
        
        while (true) {
            struct sockaddr_storage client_addr;
            socklen_t client_len = sizeof(client_addr);
            
            int client_socket = accept(local_socket, (struct sockaddr*)&client_addr, &client_len);
            if (client_socket < 0) {
                std::cerr << "Accept failed: " << strerror(errno) << std::endl;
                continue;
            }
            
            // Handle each connection in a separate thread
            std::thread([this, client_socket]() {
                this->handleConnection(client_socket);
            }).detach();
        }
        
        close(local_socket);
    }
};

void printUsage(const char* program_name) {
    std::cout << "Usage:" << std::endl;
    std::cout << "  Server mode: " << program_name << " server <listen_addr> <listen_port> <target_host> <target_port>" << std::endl;
    std::cout << "  Client mode: " << program_name << " client <local_addr> <local_port> <tunnel_host> <tunnel_port>" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  IPv4 to IPv6 tunnel server:" << std::endl;
    std::cout << "    " << program_name << " server 0.0.0.0 8080 2001:db8::1 80" << std::endl;
    std::cout << "  IPv6 to IPv4 tunnel server:" << std::endl;
    std::cout << "    " << program_name << " server :: 8080 192.168.1.100 80" << std::endl;
    std::cout << "  Tunnel client:" << std::endl;
    std::cout << "    " << program_name << " client 127.0.0.1 3128 tunnel.example.com 8080" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        printUsage(argv[0]);
        return 1;
    }
    
    std::string mode = argv[1];
    
    try {
        if (mode == "server") {
            std::string listen_addr = argv[2];
            int listen_port = std::stoi(argv[3]);
            std::string target_host = argv[4];
            int target_port = std::stoi(argv[5]);
            
            TunnelServer server(listen_addr, listen_port, target_host, target_port);
            server.start();
            
        } else if (mode == "client") {
            std::string local_host = argv[2];
            int local_port = std::stoi(argv[3]);
            std::string tunnel_host = argv[4];
            int tunnel_port = std::stoi(argv[5]);
            
            TunnelClient client(local_host, local_port, tunnel_host, tunnel_port);
            client.start();
            
        } else {
            std::cerr << "Invalid mode. Use 'server' or 'client'." << std::endl;
            printUsage(argv[0]);
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}