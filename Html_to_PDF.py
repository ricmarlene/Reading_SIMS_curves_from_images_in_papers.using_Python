import os
import pdfkit
from pathlib import Path

def html_to_pdf(source_dir, target_dir):
    """使用 pdfkit 转换 HTML 为 PDF"""
    os.makedirs(target_dir, exist_ok=True)
    
    # 配置 wkhtmltopdf 路径（需要先安装）
    config = pdfkit.configuration(wkhtmltopdf=r'C:\Program Files\wkhtmltopdf\bin\wkhtmltopdf.exe')
    
    # 添加选项：允许访问本地文件
    options = {
        'enable-local-file-access': None,
        'quiet': ''  # 静默模式，减少输出
    }
    
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            if file.lower().endswith('.html'):
                html_path = os.path.join(root, file)
                pdf_name = Path(file).with_suffix('.pdf').name
                pdf_path = os.path.join(target_dir, pdf_name)
                
                try:
                    pdfkit.from_file(
                        html_path, 
                        pdf_path, 
                        configuration=config,
                        options=options
                    )
                    print(f"转换成功: {html_path} → {pdf_path}")
                except Exception as e:
                    print(f"转换失败: {html_path} - 错误: {str(e)}")

if __name__ == "__main__":
    source_dir = input("请输入要搜索的目录: ").strip()
    while not os.path.isdir(source_dir):
        print("目录不存在！")
        source_dir = input("请重新输入: ").strip()
    
    target_dir = os.path.join(os.path.dirname(source_dir), "html_to_pdf")
    html_to_pdf(source_dir, target_dir)
