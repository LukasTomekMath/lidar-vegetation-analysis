from lxml import etree

# hlada kod pralesa a ak nie je pritomny, tak segment vymaze, zostanu iba pralesove lokality

def filter_kml_final(input_file, output_file):
    
    parser = etree.XMLParser(remove_blank_text=True)
    tree = etree.parse(input_file, parser)
    root = tree.getroot()
    
    ns = {'kml': 'http://www.opengis.net/kml/2.2'}
    
    placemarks = root.xpath('//kml:Placemark', namespaces=ns)
    
    print(f"Total Placemarks detected: {len(placemarks)}")
    
    removed_count = 0

    for placemark in placemarks:
        target_field = placemark.xpath('.//kml:SimpleData[@name="KOD_PRALES"]', namespaces=ns)
        
        should_remove = False
        
        if not target_field:
            should_remove = True
        else:
            value = target_field[0].text
            if value is None or str(value).strip() == "":
                should_remove = True
        
        if should_remove:
            parent = placemark.getparent()
            if parent is not None:
                parent.remove(placemark)
                removed_count += 1
                
    tree.write(output_file, encoding='utf-8', xml_declaration=True, pretty_print=True)
    
    print(f"Removed: {removed_count}")
    print(f"Remaining: {len(placemarks) - removed_count}")

filter_kml_final('pralesy.kml', 'filtered_output_pralesy.kml')