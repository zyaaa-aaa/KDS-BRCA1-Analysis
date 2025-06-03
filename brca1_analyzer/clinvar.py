from typing import Optional
from dataclasses import dataclass
import requests
import streamlit as st
import xml.etree.ElementTree as ET
import json
import time

@dataclass
class VariantInfo:
    """Data class for storing variant information"""
    position: int
    ref_allele: str
    alt_allele: str
    variant_name: str
    clinical_significance: str
    clinvar_id: str
    description: str
    review_status: str = ""
    condition: str = ""

class ClinVarIntegrator:
    """Integration with ClinVar database for variant classification"""
    
    def __init__(self):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.clinvar_cache = {}
        self.api_key = None  # Add your NCBI API key here for higher rate limits
    
    def search_variant(self, gene_name: str, variant: str) -> Optional[VariantInfo]:
        """
        Search for variant information in ClinVar using real API calls
        """
        cache_key = f"{gene_name}_{variant}"
        if cache_key in self.clinvar_cache:
            return self.clinvar_cache[cache_key]
        
        try:
            # First, search for the variant
            clinvar_ids = self._search_clinvar(gene_name, variant)
            
            if not clinvar_ids:
                st.warning(f"No ClinVar entries found for {gene_name} {variant}")
                return None
            
            # Fetch details for the first result
            variant_info = self._fetch_variant_details(clinvar_ids[0])
            
            if variant_info:
                self.clinvar_cache[cache_key] = variant_info
            
            return variant_info
            
        except Exception as e:
            st.error(f"ClinVar lookup failed: {str(e)}")
            return None
    
    def _search_clinvar(self, gene_name: str, variant: str) -> list:
        """
        Search ClinVar database for variant IDs
        """
        # Try different search strategies
        search_terms = [
            f"{gene_name}[gene] AND {variant}[variant name]",
            f"{gene_name} AND {variant}",
            f'"{gene_name}" AND "{variant}"'
        ]
        
        for search_term in search_terms:
            try:
                search_url = f"{self.base_url}esearch.fcgi"
                params = {
                    'db': 'clinvar',
                    'term': search_term,
                    'retmode': 'json',
                    'retmax': '10',
                    'sort': 'relevance'
                }
                
                if self.api_key:
                    params['api_key'] = self.api_key
                
                response = requests.get(search_url, params=params, timeout=10)
                response.raise_for_status()
                
                data = response.json()
                
                if 'esearchresult' in data and data['esearchresult']['idlist']:
                    return data['esearchresult']['idlist']
                
                # Add delay to respect rate limits
                time.sleep(0.34)  # ~3 requests per second
                
            except requests.exceptions.RequestException as e:
                st.warning(f"Search attempt failed for term '{search_term}': {e}")
                continue
        
        return []
    
    def _fetch_variant_details(self, clinvar_id: str) -> Optional[VariantInfo]:
        """
        Fetch detailed variant information from ClinVar
        """

        try:
            # Fetch variant details using efetch
            fetch_url = f"{self.base_url}efetch.fcgi"
            params = {
                'db': 'clinvar',
                'id': clinvar_id,
                'rettype': 'clinvarset',
                'retmode': 'xml'
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(fetch_url, params=params, timeout=10)
            response.raise_for_status()

            with open(f"clinvar_raw_{clinvar_id}.xml", "w", encoding="utf-8") as f:
                f.write(response.text)
            
            # Parse XML response
            variant_info = self._parse_clinvar_xml(response.text, clinvar_id)
            
            if variant_info:
                return variant_info
            
            # If XML parsing fails, try JSON approach
            return self._fetch_variant_summary(clinvar_id)
            
        except Exception as e:
            st.warning(f"Failed to fetch details for ClinVar ID {clinvar_id}: {e}")
            return None
    
    def _parse_clinvar_xml(self, xml_content: str, clinvar_id: str) -> Optional[VariantInfo]:
        """
        Parse XML response from ClinVar efetch - updated for actual ClinVar XML structure
        """
        try:
            root = ET.fromstring(xml_content)
            
            # Initialize variables
            variant_name = "Unknown"
            clinical_significance = "Unknown"
            description = ""
            review_status = ""
            condition = ""
            position = 0
            ref_allele = ""
            alt_allele = ""
            
            # Navigate through ClinVarSet structure
            for clinvar_set in root.findall('.//ClinVarSet'):
                
                # Get variant name from preferred name in MeasureSet
                for name_elem in clinvar_set.findall('.//MeasureSet/Measure/Name/ElementValue[@Type="preferred name"]'):
                    if name_elem.text:
                        variant_name = name_elem.text.strip()
                        break
                
                # Fallback: get variant name from any Name/ElementValue
                if variant_name == "Unknown":
                    for name_elem in clinvar_set.findall('.//MeasureSet/Measure/Name/ElementValue'):
                        if name_elem.text:
                            variant_name = name_elem.text.strip()
                            break
                
                # Get clinical significance from GermlineClassification Description
                for sig_elem in clinvar_set.findall('.//GermlineClassification/Description'):
                    if sig_elem.text:
                        clinical_significance = sig_elem.text.strip()
                        break
                
                # Fallback: get clinical significance from simple GermlineClassification tag
                if clinical_significance == "Unknown":
                    for sig_elem in clinvar_set.findall('.//GermlineClassification'):
                        if sig_elem.text and sig_elem.text.strip():
                            clinical_significance = sig_elem.text.strip()
                            break
                
                # Get review status
                for review_elem in clinvar_set.findall('.//ReviewStatus'):
                    if review_elem.text:
                        review_status = review_elem.text.strip()
                        break
                
                # Get condition from TraitSet
                for trait_elem in clinvar_set.findall('.//TraitSet/Trait/Name/ElementValue[@Type="Preferred"]'):
                    if trait_elem.text:
                        condition = trait_elem.text.strip()
                        break
                
                # Fallback: get any condition name
                if not condition:
                    for trait_elem in clinvar_set.findall('.//TraitSet/Trait/Name/ElementValue'):
                        if trait_elem.text:
                            condition = trait_elem.text.strip()
                            break
                
                # Get sequence location info
                for seq_loc in clinvar_set.findall('.//SequenceLocation'):
                    # Try to get position
                    start = seq_loc.get('start')
                    if start:
                        try:
                            position = int(start)
                        except ValueError:
                            pass
                    
                    # Get alleles if available
                    ref = seq_loc.get('referenceAllele')
                    alt = seq_loc.get('alternateAllele')
                    if ref:
                        ref_allele = ref
                    if alt:
                        alt_allele = alt
                
                # Get additional description from ObservedData
                for desc_elem in clinvar_set.findall('.//ObservedData/Attribute[@Type="Description"]'):
                    if desc_elem.text and len(desc_elem.text.strip()) > len(description):
                        description = desc_elem.text.strip()
                
                break  # Use first ClinVarSet
            
            # Build final description if not found in ObservedData
            if not description:
                desc_parts = []
                if clinical_significance != "Unknown":
                    desc_parts.append(f"{clinical_significance} variant")
                if condition:
                    desc_parts.append(f"associated with {condition}")
                if review_status:
                    desc_parts.append(f"(Review status: {review_status})")
                
                description = " ".join(desc_parts) if desc_parts else "Variant found in ClinVar"
            
            # Clean up clinical significance - capitalize first letter
            if clinical_significance != "Unknown":
                clinical_significance = clinical_significance.capitalize()
            
            return VariantInfo(
                position=position,
                ref_allele=ref_allele,
                alt_allele=alt_allele,
                variant_name=variant_name,
                clinical_significance=clinical_significance,
                clinvar_id=clinvar_id,
                description=description,
                review_status=review_status,
                condition=condition
            )
            
        except ET.ParseError as e:
            st.warning(f"XML parsing error: {e}")
            return None
        except Exception as e:
            st.warning(f"Error parsing ClinVar XML: {e}")
            return None
    
    def _fetch_variant_summary(self, clinvar_id: str) -> Optional[VariantInfo]:
        """
        Fetch variant summary using esummary as fallback
        """
        try:
            summary_url = f"{self.base_url}esummary.fcgi"
            params = {
                'db': 'clinvar',
                'id': clinvar_id,
                'retmode': 'json'
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(summary_url, params=params, timeout=10)
            response.raise_for_status()
            
            data = response.json()
            
            if 'result' in data and clinvar_id in data['result']:
                result = data['result'][clinvar_id]
                
                return VariantInfo(
                    position=0,
                    ref_allele="",
                    alt_allele="",
                    variant_name=result.get('title', 'Unknown variant'),
                    clinical_significance=result.get('clinical_significance', 'Unknown'),
                    clinvar_id=clinvar_id,
                    description=result.get('summary', 'No description available'),
                    review_status=result.get('review_status', ''),
                    condition=result.get('condition', '')
                )
            
            return None
            
        except Exception as e:
            st.warning(f"Failed to fetch summary for ClinVar ID {clinvar_id}: {e}")
            return None
    
    def search_multiple_variants(self, gene_name: str, variants: list) -> dict:
        """
        Search for multiple variants efficiently
        """
        results = {}
        
        for variant in variants:
            try:
                result = self.search_variant(gene_name, variant)
                results[variant] = result
                
                # Add delay between requests to respect rate limits
                time.sleep(0.34)
                
            except Exception as e:
                st.warning(f"Failed to search variant {variant}: {e}")
                results[variant] = None
        
        return results
    
    def get_gene_variants(self, gene_name: str, max_results: int = 50) -> list:
        """
        Get all variants for a specific gene
        """
        try:
            search_url = f"{self.base_url}esearch.fcgi"
            params = {
                'db': 'clinvar',
                'term': f"{gene_name}[gene]",
                'retmode': 'json',
                'retmax': str(max_results),
                'sort': 'relevance'
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(search_url, params=params, timeout=10)
            response.raise_for_status()
            
            data = response.json()
            
            if 'esearchresult' in data and data['esearchresult']['idlist']:
                clinvar_ids = data['esearchresult']['idlist']
                
                variants = []
                for clinvar_id in clinvar_ids[:10]:  # Limit to first 10 for performance
                    variant_info = self._fetch_variant_details(clinvar_id)
                    if variant_info:
                        variants.append(variant_info)
                    
                    time.sleep(0.34)  # Rate limiting
                
                return variants
            
            return []
            
        except Exception as e:
            st.error(f"Failed to get variants for gene {gene_name}: {e}")
            return []