import heapq
from collections import defaultdict

class TransportGraph:
    def __init__(self):
        self.stations = []
        self.station_index = {}
        self.adjacency_list = defaultdict(list)
        self.edges = []
        self.distance_matrix = None
        
    def add_station(self, name):
        if name not in self.station_index:
            self.station_index[name] = len(self.stations)
            self.stations.append(name)
            
    def add_route(self, from_station, to_station, weight):
        self.add_station(from_station)
        self.add_station(to_station)
        
        from_idx = self.station_index[from_station]
        to_idx = self.station_index[to_station]
        
        # Добавляем ребро в список смежности (для Дейкстры и Беллмана-Форда)
        self.adjacency_list[from_idx].append((to_idx, weight))
        self.adjacency_list[to_idx].append((from_idx, weight))
        
        # Добавляем ребро в список всех рёбер (для Крускала)
        self.edges.append((from_idx, to_idx, weight))
        
    def shortest_path_Dijkstra(self, start):
        """Алгоритм Дейкстры для поиска кратчайших путей без отрицательных весов"""
        if start not in self.station_index:
            return {}
            
        start_idx = self.station_index[start]
        n = len(self.stations)
        distances = [float('inf')] * n
        distances[start_idx] = 0
        
        pq = [(0, start_idx)]
        
        while pq:
            current_dist, current_node = heapq.heappop(pq)
            
            if current_dist > distances[current_node]:
                continue
                
            for neighbor, weight in self.adjacency_list[current_node]:
                distance = current_dist + weight
                if distance < distances[neighbor]:
                    distances[neighbor] = distance
                    heapq.heappush(pq, (distance, neighbor))
        
        # Преобразуем в словарь с именами станций
        result = {}
        for i, dist in enumerate(distances):
            result[self.stations[i]] = dist
            
        return result
    
    def shortest_path_BellmanFord(self, start):
        """Алгоритм Беллмана-Форда для работы с отрицательными весами"""
        if start not in self.station_index:
            return {}
            
        start_idx = self.station_index[start]
        n = len(self.stations)
        distances = [float('inf')] * n
        distances[start_idx] = 0
        
        # Создаем список всех направленных рёбер
        edges = []
        for u in range(n):
            for v, weight in self.adjacency_list[u]:
                edges.append((u, v, weight))
        
        # Релаксация рёбер n-1 раз
        for _ in range(n - 1):
            updated = False
            for u, v, weight in edges:
                if distances[u] != float('inf') and distances[u] + weight < distances[v]:
                    distances[v] = distances[u] + weight
                    updated = True
            if not updated:
                break
        
        # Проверка на отрицательные циклы
        for u, v, weight in edges:
            if distances[u] != float('inf') and distances[u] + weight < distances[v]:
                raise ValueError("Граф содержит отрицательный цикл")
        
        # Преобразуем в словарь с именами станций
        result = {}
        for i, dist in enumerate(distances):
            result[self.stations[i]] = dist
            
        return result
    
    def all_pairs_shortest_paths_FloydWarshall(self):
        """Алгоритм Флойда-Уоршелла для всех пар вершин"""
        n = len(self.stations)
        
        # Инициализация матрицы расстояний
        dist = [[float('inf')] * n for _ in range(n)]
        
        for i in range(n):
            dist[i][i] = 0
            
        for u in range(n):
            for v, weight in self.adjacency_list[u]:
                dist[u][v] = weight
        
        # Алгоритм Флойда-Уоршелла
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    if dist[i][k] != float('inf') and dist[k][j] != float('inf'):
                        if dist[i][j] > dist[i][k] + dist[k][j]:
                            dist[i][j] = dist[i][k] + dist[k][j]
        
        self.distance_matrix = dist
        return dist

class NetworkOptimizer:
    def __init__(self, graph):
        self.graph = graph
        
    def minimum_spanning_tree_Kruskal(self):
        """Алгоритм Крускала для построения минимального остовного дерева"""
        n = len(self.graph.stations)
        
        # Сортируем рёбра по весу
        edges = sorted(self.graph.edges, key=lambda x: x[2])
        
        parent = list(range(n))
        rank = [0] * n
        
        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]
        
        def union(x, y):
            root_x, root_y = find(x), find(y)
            if root_x == root_y:
                return False
            if rank[root_x] < rank[root_y]:
                parent[root_x] = root_y
            elif rank[root_x] > rank[root_y]:
                parent[root_y] = root_x
            else:
                parent[root_y] = root_x
                rank[root_x] += 1
            return True
        
        mst_edges = []
        total_weight = 0
        
        for u, v, weight in edges:
            if union(u, v):
                mst_edges.append((self.graph.stations[u], self.graph.stations[v], weight))
                total_weight += weight
                if len(mst_edges) == n - 1:
                    break
        
        return mst_edges, total_weight
    
    def minimum_spanning_tree_Prim(self):
        """Алгоритм Прима для построения минимального остовного дерева"""
        n = len(self.graph.stations)
        visited = [False] * n
        mst_edges = []
        total_weight = 0
        
        # Начинаем с вершины 0
        pq = []
        heapq.heappush(pq, (0, 0, -1))  # (weight, current, parent)
        
        while pq and len(mst_edges) < n - 1:
            weight, current, parent = heapq.heappop(pq)
            
            if visited[current]:
                continue
                
            visited[current] = True
            
            if parent != -1:
                mst_edges.append((self.graph.stations[parent], self.graph.stations[current], weight))
                total_weight += weight
            
            for neighbor, edge_weight in self.graph.adjacency_list[current]:
                if not visited[neighbor]:
                    heapq.heappush(pq, (edge_weight, neighbor, current))
        
        return mst_edges, total_weight

def create_city_network():
    """Создание тестовой транспортной сети города"""
    graph = TransportGraph()
    
    # Добавляем станции
    stations = [
        "Центр. вокзал", "Пл. Ленина", "Автостанция-1", "Университет",
        "Вокзал-Северный", "Автостанция-2", "Парковая", "Портовая"
    ]
    
    # Добавляем маршруты согласно матрице
    routes = [
        ("Центр. вокзал", "Пл. Ленина", 1.2),
        ("Центр. вокзал", "Автостанция-1", 3.5),
        ("Центр. вокзал", "Вокзал-Северный", 4.0),
        
        ("Пл. Ленина", "Автостанция-1", 2.1),
        ("Пл. Ленина", "Университет", 0.8),
        ("Пл. Ленина", "Автостанция-2", 3.0),
        ("Пл. Ленина", "Парковая", 1.5),
        
        ("Автостанция-1", "Автостанция-2", 5.4),
        ("Автостанция-1", "Портовая", 4.2),
        
        ("Университет", "Вокзал-Северный", 2.7),
        ("Университет", "Автостанция-2", 1.9),
        
        ("Вокзал-Северный", "Парковая", 3.3),
        ("Вокзал-Северный", "Портовая", 6.0),
        
        ("Автостанция-2", "Парковая", 2.5),
    ]
    
    for station in stations:
        graph.add_station(station)
    
    for from_station, to_station, weight in routes:
        graph.add_route(from_station, to_station, weight)
    
    return graph

def print_matrix(matrix, stations):
    """Красивая печать матрицы расстояний"""
    print("\nМатрица кратчайших расстояний:")
    print(" " * 15, end="")
    for station in stations:
        print(f"{station[:8]:>8}", end="")
    print()
    
    for i, station in enumerate(stations):
        print(f"{station:15}", end="")
        for j in range(len(stations)):
            if matrix[i][j] == float('inf'):
                print(f"{'-':>8}", end="")
            else:
                print(f"{matrix[i][j]:>8.1f}", end="")
        print()

def main():
    print("=" * 60)
    print("ОПТИМИЗАЦИЯ ГОРОДСКОЙ ТРАНСПОРТНОЙ СЕТИ")
    print("=" * 60)
    
    # Создаем транспортную сеть города
    city_graph = create_city_network()
    optimizer = NetworkOptimizer(city_graph)
    
    # 1. Поиск кратчайших путей от Центрального вокзала (Дейкстра)
    print("\n1. КРАТЧАЙШИЕ ПУТИ ОТ ЦЕНТРАЛЬНОГО ВОКЗАЛА (Дейкстра):")
    dijkstra_result = city_graph.shortest_path_Dijkstra("Центр. вокзал")
    for station, distance in sorted(dijkstra_result.items()):
        print(f"   {station:20} → {distance:5.1f} км")
    
    # 2. Проверка на отрицательные циклы (Беллман-Форд)
    print("\n2. ПРОВЕРКА ОТРИЦАТЕЛЬНЫХ ЦИКЛОВ (Беллман-Форд):")
    try:
        bellman_result = city_graph.shortest_path_BellmanFord("Центр. вокзал")
        print("   ✓ Отрицательных циклов не обнаружено")
    except ValueError as e:
        print(f"   ✗ {e}")
    
    # 3. Матрица всех кратчайших путей (Флойд-Уоршелл)
    floyd_matrix = city_graph.all_pairs_shortest_paths_FloydWarshall()
    print_matrix(floyd_matrix, city_graph.stations)
    
    # 4. Минимальная остовное дерево (Крускал)
    print("\n4. МИНИМАЛЬНАЯ ТРАНСПОРТНАЯ СЕТЬ (Крускал):")
    mst_kruskal, total_kruskal = optimizer.minimum_spanning_tree_Kruskal()
    print("   Оптимальные соединения:")
    for u, v, weight in mst_kruskal:
        print(f"   {u:15} ↔ {v:15} ({weight:4.1f} км)")
    print(f"   Суммарная длина: {total_kruskal:.1f} км")
    
    # 5. Минимальная остовное дерево (Прим)
    print("\n5. МИНИМАЛЬНАЯ ТРАНСПОРТНАЯ СЕТЬ (Прим):")
    mst_prim, total_prim = optimizer.minimum_spanning_tree_Prim()
    print("   Оптимальные соединения:")
    for u, v, weight in mst_prim:
        print(f"   {u:15} ↔ {v:15} ({weight:4.1f} км)")
    print(f"   Суммарная длина: {total_prim:.1f} км")
    
    # Анализ результатов
    print("\n" + "=" * 60)
    print("АНАЛИЗ РЕЗУЛЬТАТОВ:")
    print("=" * 60)
    
    # Самые удаленные станции от центра
    farthest = sorted(dijkstra_result.items(), key=lambda x: x[1], reverse=True)[:3]
    print(f"Самые удаленные станции от Центрального вокзала:")
    for station, dist in farthest:
        print(f"   • {station}: {dist:.1f} км")
    
    # Критические соединения в минимальной сети
    print(f"\nКритические соединения минимальной сети:")
    for u, v, weight in mst_kruskal:
        if weight > 3.0:
            print(f"   • {u} ↔ {v} ({weight:.1f} км) - важное дальнее соединение")

if __name__ == "__main__":

    main()
