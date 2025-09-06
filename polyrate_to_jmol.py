import re
import os

class PolyrateJmolConverter:
    """
    一个用于解析 Polyrate 输出文件，提取优化几何结构、振动频率和简正模式，
    并将其转换为 Jmol 可识别的 Molden 格式文件的工具。
    (版本 14：修复原子信息解析问题)
    """

    SPECIES_CONFIG = {
        'reactant1':      {'calc_header': 'Reactant 1',       'param_header': 'Reactant #1 Parameters:'},
        'reactant2':      {'calc_header': 'Reactant 2',       'param_header': 'Reactant #2 Parameters:'},
        'product1':       {'calc_header': 'Product 1',        'param_header': 'Product #1 Parameters:'},
        'product2':       {'calc_header': 'Product 2',        'param_header': 'Product #2 Parameters:'},
        'reactants_well': {'calc_header': 'Reactants well',   'param_header': 'WELLR Parameters:'},
        'products_well':  {'calc_header': 'Products well',    'param_header': 'WELLP Parameters:'},
        'saddle_point':   {'calc_header': 'Saddle point',     'param_header': 'Starting point Parameters:'}
    }

    # Angstrom 转 Bohr 的转换因子
    ANG_TO_BOHR = 1.8897259886

    def __init__(self, input_filepath):
        self.input_filepath = input_filepath
        try:
            with open(input_filepath, 'r', encoding='utf-8') as f:
                self.content = f.read()
            print(f"成功读取文件: {input_filepath}")
        except FileNotFoundError:
            print(f"错误: 文件未找到 '{input_filepath}'")
            self.content = None
            return

        self.atom_info = self._parse_atom_info()  # 现在解析原子符号和质量
        self.calculation_blocks = self._split_into_calculation_blocks()

    def _parse_atom_info(self):
        atom_info = {}
        pattern = re.compile(r"Atomic information:.*?\n\n(.*?)\n\n---", re.DOTALL)
        match = pattern.search(self.content)
        if not match:
            print("警告: 未在文件中找到 'Atomic information' 部分。")
            return {}
            
        # 提取原子信息表格内容
        atom_table = match.group(1).strip()
        lines = atom_table.split('\n')
        
        # 跳过表头行（如果有）
        start_index = 0
        for i, line in enumerate(lines):
            if line.strip().startswith("Atom") or "Atomic number" in line:
                start_index = i + 1
                break
        
        # 解析原子数据行
        for line in lines[start_index:]:
            parts = re.split(r'\s{2,}', line.strip())  # 基于多个空格分割
            if len(parts) >= 4 and parts[0].isdigit():
                atom_id = parts[0]
                # 原子符号可能在第三或第四列
                symbol = parts[2] if not parts[2].replace('.', '', 1).isdigit() else parts[3]
                
                # 质量可能在最后一列
                mass_str = parts[-1].strip()
                try:
                    mass = float(mass_str)
                    atom_info[atom_id] = {'symbol': symbol, 'mass': mass}
                    print(f"解析原子: ID={atom_id}, 符号={symbol}, 质量={mass:.6f} amu")
                except ValueError:
                    print(f"警告: 无法解析原子质量 '{mass_str}'，跳过原子 {atom_id}")
        
        print(f"成功解析到 {len(atom_info)} 个原子的信息（包括质量）。")
        return atom_info

    def _split_into_calculation_blocks(self):
        blocks = {}
        delimiter_pattern = r'(\*{20,}.*?\*{20,})'
        chunks = re.split(delimiter_pattern, self.content)
        
        for i in range(1, len(chunks), 2):
            header = chunks[i].replace('*', '').strip()
            content_block = chunks[i+1]
            blocks[header] = content_block
        
        print(f"文件被成功分割成 {len(blocks)} 个计算块。")
        return blocks

    def _parse_final_geometry(self, calc_block):
        pattern = re.compile(
            r"Final geometry in unscaled Cartesians \(angstroms\)\n\n Atom\s+X\s+Y\s+Z\n\n(.*?)(?:\n\n Final derivatives|\n\n  V =)",
            re.DOTALL
        )
        match = pattern.search(calc_block)
        if not match:
            print("警告: 未找到最终几何结构。")
            return None
        
        geometry = {}
        lines = match.group(1).strip().split('\n')
        for line in lines:
            parts = line.split()
            if len(parts) == 4 and parts[0].isdigit():
                atom_id, x, y, z = parts
                geometry[atom_id] = [float(x), float(y), float(z)]
        print(f"成功解析到最终几何结构，包含 {len(geometry)} 个原子。")
        return geometry
        
    def _parse_initial_geometry(self, param_header_text):
        pattern = re.compile(
            rf"{re.escape(param_header_text)}.*?Geom.*?X\s+Y\s+Z\n(.*?)(?=\n\n\s+\*\* All coordinates)",
            re.DOTALL
        )
        match = pattern.search(self.content)
        if not match:
            print(f"警告: 未在参数部分 '{param_header_text}' 中找到初始几何结构。")
            return None
            
        geometry = {}
        lines = match.group(1).strip().split('\n')
        for line in lines:
            parts = line.split()
            if len(parts) == 5 and parts[0].isdigit():
                atom_id, _, x, y, z = parts
                geometry[atom_id] = [float(x), float(y), float(z)]
        print(f"成功从参数部分解析到初始几何结构，包含 {len(geometry)} 个原子。")
        return geometry

    def _parse_frequencies_and_modes(self, calc_block, num_atoms):
        """
        从计算块中解析振动频率和对应的简正模式向量。
        """
        pattern = re.compile(
            r"Frequencies and normalized eigenvector components(.*?)(?=\n -{5,}|\Z)", 
            re.DOTALL
        )
        match = pattern.search(calc_block)
        if not match:
            print("警告: 未找到 'Frequencies and normalized eigenvector components' 部分")
            return []

        print("找到振动数据部分")
        vibrations = []
        lines = match.group(1).strip().split('\n')
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("Mode"):
                print(f"\n处理模式块: {line}")
                i += 1
                if i >= len(lines):
                    print("警告: 模式块不完整，缺少频率行")
                    break
                
                # 解析频率行
                freq_line = lines[i].strip()
                print(f"频率行: {freq_line}")
                if freq_line.startswith("omega"):
                    freq_parts = [p for p in freq_line.split() if p]
                    freq_values = []
                    j = 2  # 跳过 "omega (cm**-1)"
                    while j < len(freq_parts):
                        try:
                            if freq_parts[j].replace('.', '', 1).replace('-', '', 1).isdigit():
                                if j + 1 < len(freq_parts) and freq_parts[j + 1] == 'i':
                                    freq_values.append(-float(freq_parts[j]))  # 虚频
                                    j += 2
                                else:
                                    freq_values.append(float(freq_parts[j]))
                                    j += 1
                            else:
                                j += 1
                        except ValueError as e:
                            print(f"警告: 无法解析频率值 '{freq_parts[j]}'，错误: {e}")
                            j += 1
                    print(f"解析到的频率: {freq_values}")
                    
                    if not freq_values:
                        print("警告: 未从频率行中提取到有效频率")
                        i += 1
                        continue

                    # 跳过频率行后的无关行，寻找向量数据
                    i += 1
                    vector_sets = [{} for _ in freq_values]
                    atom_count = 0
                    
                    while i < len(lines) and atom_count < num_atoms:
                        vector_line = lines[i].strip()
                        print(f"检查行: {vector_line}")
                        parts = vector_line.split()
                        # 假设向量行以数字开头，后跟坐标数据
                        if parts and parts[0].isdigit() and len(parts) >= 1 + len(freq_values) * 3:
                            atom_id = parts[0]
                            print(f"解析向量行: {vector_line}")
                            try:
                                for k in range(len(freq_values)):
                                    start_idx = 1 + k * 3
                                    vector_sets[k][atom_id] = [float(parts[start_idx + m]) for m in range(3)]
                                atom_count += 1
                            except (ValueError, IndexError) as e:
                                print(f"警告: 无法解析向量行 '{vector_line}'，错误: {e}")
                        i += 1
                    
                    if atom_count < num_atoms:
                        print(f"警告: 仅解析到 {atom_count} 个原子的向量，预期 {num_atoms} 个")
                    
                    # 将频率和向量组合
                    for freq, vectors in zip(freq_values, vector_sets):
                        if vectors and len(vectors) == num_atoms:
                            vibrations.append({'freq': freq, 'vectors': vectors})
                        else:
                            print(f"警告: 频率 {freq} 的向量数据不完整，仅包含 {len(vectors)} 个原子")

            else:
                i += 1

        print(f"\n总共解析到 {len(vibrations)} 个振动模式")
        return vibrations

    def _convert_mass_weighted_to_cartesian(self, vibrations):
        """
        将质量加权的简正模式转换为笛卡尔位移向量并进行归一化
        """
        if not vibrations or not self.atom_info:
            return vibrations
            
        print("\n开始转换质量加权坐标为笛卡尔位移...")
        
        for vibration in vibrations:
            vectors = vibration['vectors']
            total_square_sum = 0.0
            
            # 第一步：质量反加权
            for atom_id, disp in vectors.items():
                if atom_id in self.atom_info:
                    mass = self.atom_info[atom_id]['mass']
                    sqrt_mass = mass**0.5
                    # 质量反加权: Δr = q / √m
                    cart_disp = [d / sqrt_mass for d in disp]
                    vectors[atom_id] = cart_disp
                    
                    # 累加平方和用于归一化
                    for d_val in cart_disp:
                        total_square_sum += d_val**2
                else:
                    print(f"警告: 原子ID {atom_id} 的质量信息缺失，跳过质量反加权")
            
            # 第二步：归一化 (确保 ∑(Δr)^2 = 1)
            norm_factor = total_square_sum**0.5
            if norm_factor < 1e-10:
                print(f"警告: 振动模式 {vibration['freq']} 的向量模长接近零，跳过归一化")
                continue
                
            print(f"模式 {vibration['freq']:.2f} cm⁻¹ 的归一化因子: {norm_factor:.6f}")
            
            for atom_id in vectors:
                vectors[atom_id] = [d / norm_factor for d in vectors[atom_id]]
                
        return vibrations

    def _write_molden_file(self, output_filename, title, geometry, vibrations):
        """
        生成符合Molden格式的文件，Jmol可以直接加载并可视化振动模式
        """
        if not geometry:
            print(f"错误: 无法为 '{title}' 生成 Molden 文件，因为缺少几何信息。")
            return

        # 过滤掉频率绝对值小于10的振动模式
        valid_vibrations = [v for v in vibrations if abs(v['freq']) > 10.0]
        
        with open(output_filename, 'w', encoding='utf-8') as f:
            # 写入文件头
            f.write("[Molden Format]\n")
            f.write(f"# Generated from Polyrate output for {title}\n")
            f.write("\n")

            # 写入几何结构部分 (使用[GEOMETRIES] XYZ格式，单位Å)
            f.write("[GEOMETRIES] XYZ\n")
            sorted_atom_ids = sorted(geometry.keys(), key=int)
            num_atoms = len(sorted_atom_ids)
            
            # XYZ格式第一部分：原子数
            f.write(f"{num_atoms}\n")
            # XYZ格式第二部分：标题行（使用物种名称）
            f.write(f"{title}\n")
            # XYZ格式第三部分：原子坐标（单位Å）
            for atom_id in sorted_atom_ids:
                if atom_id in self.atom_info:
                    symbol = self.atom_info[atom_id]['symbol']
                else:
                    symbol = 'X'
                    print(f"警告: 原子ID {atom_id} 的符号信息缺失，使用 'X' 代替")
                
                x, y, z = geometry[atom_id]  # 注意这里已经是Å单位
                f.write(f"{symbol:>2s}  {x:>15.8f} {y:>15.8f} {z:>15.8f}\n")            
            # 写入频率部分
            if valid_vibrations:
                f.write("\n[FREQ]\n")
                for iv,vib in enumerate(valid_vibrations):
                    f.write(f"  {vib['freq']:>10.4f}")
                    if iv != (len(valid_vibrations)-1):
                        f.write(f"\n")
            
            # 写入参考坐标（与原子坐标相同，转换为Bohr）
            f.write("\n[FR-COORD]\n")
            for atom_id in sorted_atom_ids:
                if atom_id in self.atom_info:
                    symbol = self.atom_info[atom_id]['symbol']
                else:
                    symbol = 'X'
                
                x, y, z = geometry[atom_id]
                # 转换为Bohr单位
                x_bohr = x * self.ANG_TO_BOHR
                y_bohr = y * self.ANG_TO_BOHR
                z_bohr = z * self.ANG_TO_BOHR
                f.write(f"  {symbol:>2s}  {x_bohr:>15.8f} {y_bohr:>15.8f} {z_bohr:>15.8f}\n")
            
            # 写入振动向量部分 (已经是笛卡尔位移，单位bohr)
            if valid_vibrations:
                f.write("\n[FR-NORM-COORD]\n")
                for vib in valid_vibrations:
                    f.write(f" vibration\t{vib['freq']:.2f} cm-1\n")
                    for atom_id in sorted_atom_ids:
                        dx, dy, dz = vib['vectors'].get(atom_id, [0.0, 0.0, 0.0])
                        f.write(f"  {dx:>15.8f} {dy:>15.8f} {dz:>15.8f}\n")
        
        print(f"成功创建 Molden 文件: {output_filename} (包含 {len(valid_vibrations)} 个振动模式)")

    def process_all(self):
        if not self.content:
            print("因文件读取失败，处理中止。")
            return
        
        for name, config in self.SPECIES_CONFIG.items():
            print(f"\n--- 正在处理: {config['calc_header']} ---")
            
            calc_block = self.calculation_blocks.get(config['calc_header'])
            if not calc_block:
                print(f"警告: 未找到物种 '{config['calc_header']}' 的计算块。・・・")
                continue

            geometry = self._parse_final_geometry(calc_block)
            
            if not geometry:
                print(f"未找到 '{config['calc_header']}' 的最终几何结构，尝试从参数部分解析初始几何...")
                geometry = self._parse_initial_geometry(config['param_header'])
            
            if not geometry:
                print(f"错误: 无法为 '{config['calc_header']}' 找到任何几何结构。跳过。")
                continue
            
            num_atoms = len(geometry)
            vibrations = self._parse_frequencies_and_modes(calc_block, num_atoms)
            
            # 关键修改：将质量加权坐标转换为笛卡尔位移并进行归一化
            vibrations = self._convert_mass_weighted_to_cartesian(vibrations)
            
            output_filename = f"{name}_vibrations.molden"
            self._write_molden_file(output_filename, config['calc_header'], geometry, vibrations)

if __name__ == "__main__":
    input_file = "test-ch5o2.fu6"
    if not os.path.exists(input_file):
        print(f"错误: 输入文件 '{input_file}' 不在当前目录中。")
    else:
        converter = PolyrateJmolConverter(input_file)
        converter.process_all()
